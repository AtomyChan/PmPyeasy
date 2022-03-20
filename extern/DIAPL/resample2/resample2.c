/*========================================================*/
/*                                                        */
/*  resample2.c         version 1.5.1   2005.05.12        */
/*                                                        */
/*  Original source code:                                 */
/*  Copyright (C) 2005 by Przemek Wozniak                 */
/*  wozniak@lanl.gov                                      */
/*                                                        */
/*  Change bad column to background with noise            */
/*  by Igor Soszynski & Karol Zebrun                      */
/*                                                        */
/*  Modifications:                                        */
/*  Copyright (C) 2005 by Wojtek Pych, CAMK PAN           */
/*  pych@camk.edu.pl                                      */
/*                                                        */
/*  Written for GNU project C and C++ Compiler            */
/*                                                        */
/* Using coefficients of the coordinate transformation    */
/* stored in a binary file interpolates a given image     */
/* onto the grid of a reference image.                    */
/* Bicubic spline interpolation is used.                  */
/*                                                        */
/*========================================================*/

/***************************************************************************/
/*                                                                         */
/*   This program is free software; you can redistribute it and/or modify  */
/*   it under the terms of the GNU General Public License as published by  */
/*   the Free Software Foundation; either version 2 of the License, or     */
/*   (at your option) any later version.                                   */
/*                                                                         */
/*   This program is distributed in the hope that it will be useful,       */
/*   but WITHOUT ANY WARRANTY; without even the implied warranty of        */
/*   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         */
/*   GNU General Public License for more details.                          */
/*                                                                         */
/*   You should have received a copy of the GNU General Public License     */
/*   along with this program; if not, write to the                         */
/*   Free Software Foundation, Inc.,                                       */
/*   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             */
/*                                                                         */
/***************************************************************************/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "errmess.h"
#include "pfitsio1.h"

#include "defs.h"

#include "bicspl.h"
#include "fix_cosmics.h"
#include "hedit.h"
#include "rm_bad_col.h"
#include "res_bad_col.h"

/*--------------------------------------------------------*/
void read_params(char *name, char *instrname, RINGPAR *par)
{
        char    line[201],
                key[201],
                val[201];
        int     i;
        FILE    *inf;

  if (!(inf=fopen(name, "r"))) errmess(name);

    for (i=0; !feof(inf); i++)
  {
    fgets(line, 200, inf);
    if (feof(inf)) break;

    sscanf(line, "%s = %s", key, val);

    if (!strcasecmp(key, "END")) break;
    else if (!strcasecmp(key, "FIX_COSMICS")) sscanf(val, "%hd", &par->fixcosmics);
    else if (!strcasecmp(key, "MAX_GRAD"))    sscanf(val, "%f", &par->maxgrad);
    else if (!strcasecmp(key, "N_SIG_COSMICS")) sscanf(val, "%f", &par->nsigcosmics);
    else if (!strcasecmp(key, "FIX_RINGS"))   sscanf(val, "%d", &par->fix);
    else if (!strcasecmp(key, "N_SIG_RINGS")) sscanf(val, "%f", &par->nsig);
    else if (!strcasecmp(key, "N_SIG_RM"))    sscanf(val, "%f", &par->nsigrm);
    else if (!strcasecmp(key, "GROW_RAD"))    sscanf(val, "%d", &par->growrad);
    else if (!strcasecmp(key, "VERBOSE"))     sscanf(val, "%d", &par->verbose);
    else printf("WARNING: Unknown parameter %s\n", key);
  }

  fclose(inf);

  if (par->verbose > 1)
  {
    printf("Reading %s\n", name);
    printf("----------\n");
    printf("FIX_COSMICS =   %hd\n", par->fixcosmics);
    printf("MAX_GRAD =      %g\n", par->maxgrad);
    printf("N_SIG_COSMICS = %g\n", par->nsigcosmics);
    printf("FIX_RINGS =     %d\n", par->fix);
    printf("N_SIG_RINGS =   %g\n", par->nsig);
    printf("N_SIG_RM =      %g\n", par->nsigrm);
    printf("GROW_RAD =      %d\n", par->growrad);
    printf("VERBOSE =       %d\n", par->verbose);
    printf("----------\n");
  }

  if (!(inf=fopen(instrname, "r"))) errmess(instrname);

  for (i=0; !feof(inf); i++)
  {
    fgets(line, 200, inf);
    if (feof(inf)) break;

    sscanf(line, "%s = %s", key, val);

    if (!strcasecmp(key, "END")) break;
    else if (!strcasecmp(key, "SAT_LEVEL")) sscanf(val, "%f", &par->sat);
    else if (!strcasecmp(key, "MIN_LEVEL")) sscanf(val, "%f", &par->min);
    else if (!strcasecmp(key, "GAIN"))      sscanf(val, "%f", &par->gain);
    else if (!strcasecmp(key, "NX")) continue;
    else if (!strcasecmp(key, "NY")) continue;
    else if (!strcasecmp(key, "MARG")) continue;
    else if (!strcasecmp(key, "FITS")) continue;
    else printf("WARNING: Unknown parameter %s\n", key);
  }

  fclose(inf);

  if (par->verbose > 1)
  {
    printf("Reading %s\n", instrname);
    printf("-----------\n");
    printf("SAT_LEVEL = %g\n", par->sat);
    printf("MIN_LEVEL = %g\n", par->min);
    printf("GAIN =      %g\n", par->gain);
    printf("-----------\n");
  }

  return;
}
/*--------------------------------------------------------*/
int readcoeff(char *coefname, double **coefx, double **coefy)
{
        int     nterm;
        double  *cx, *cy;
        FILE    *inpf;

  if (!(inpf=fopen(coefname, "r"))) errmess(coefname);

  fread(&nterm, sizeof(int), 1, inpf);

  if (!(cx=(double *)calloc(nterm, sizeof(double)))) errmess("calloc(cx)");
  if (!(cy=(double *)calloc(nterm, sizeof(double)))) errmess("calloc(cy)");

  fread(cx, sizeof(double), nterm, inpf);
  fread(cy, sizeof(double), nterm, inpf);

  fclose(inpf);

  *coefx=cx;
  *coefy=cy;

  return(nterm);
}
/*--------------------------------------------------------*/
int main(int argc, char *argv[])
{
        char    **header,
                *parfname,
                *instrname,
                *im1name, *im2name,
                *coefname;
        int     nx, ny, nterm, ndeg, npix, i, *col, hsize;
        float   *im1, *im2;
        double  *coeffx, *coeffy;
        RINGPAR par;

/* IO stuff */
  if (argc != 6)
  {
    printf("\n\tUSAGE: resample2  param_file instrument_file coeff_file ");
    printf("input_image output_image\n\n");
    exit(-1);
  }

  parfname = argv[1];
  instrname= argv[2];
  coefname = argv[3];
  im1name  = argv[4];
  im2name  = argv[5];

/* read in parameters */
  read_params(parfname, instrname, &par);

/* read in coefficients */
  nterm=readcoeff(coefname, &coeffx, &coeffy);
  ndeg = (int)floor(-1.5+sqrt(0.25+2.0*nterm));
  if (par.verbose)
  {
    printf("%s\n", coefname);
    printf("nterm= %d   ndeg= %d\n", nterm, ndeg);
    if (par.verbose > 1)
      for (i=0; i<nterm; i++)
        printf("coeffx[%d]= %g   coeffy[%d]= %g\n", i, coeffx[i], i, coeffy[i]);
  }

/* read in the image */
  im1=read_FITS_2D1file(im1name, 's', &hsize, &header, &nx, &ny);
  npix=nx*ny;

  if (par.fixcosmics)
    fix_cosmics(im1, nx, ny, par.maxgrad, par.nsigcosmics, par.gain);

/* change bad column to background with noise (Igor Soszynski & Karol Zebrun) */
  if (!(col=(int *)malloc(nx*sizeof(int)))) errmess("malloc(col)");

  rm_bad_col(par, im1, nx, ny, col);

  if (!(im2=(float *)malloc(npix*sizeof(float)))) errmess("malloc(im2)");

/* the main procedure doing bicubic spline interpolations for all points */
  bicspl(im1, im2, nx, ny, coeffx, coeffy, ndeg, par);

  free(im1);

/* Restores bad column (Igor Soszynski & Karol Zebrun) */
  res_bad_col(par, im2, nx, ny, col, coeffx, ndeg);

  free(col);

/* change pixel type to float */
  hedit(hsize, header, "BITPIX  ", 'i', "Bits per pixel",
                          "                 %3d", -32);
  hsize=del_header_card(hsize, &header, "BZERO");
  hsize=del_header_card(hsize, &header, "BSCALE");

/* write interpolated image with new header */
  write_FITS_2D1file(im2name, hsize, header, npix, 4, im2);

  free(coeffx);
  free(coeffy);

  for (i=0; i<hsize; i++) free(header[i]);
  free(header);
  free(im2);

  return(0);
}
/*** END ***/
