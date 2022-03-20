/*========================================================*/
/*                                                        */
/*  sfind.c             version 1.6.3   2006.10.23        */
/*                                                        */
/*  Original source code:                                 */
/*  Copyright (C) 2006 by Przemek Wozniak                 */
/*  wozniak@lanl.gov                                      */
/*                                                        */
/*  Modifications:                                        */
/*  Copyright (C) 2006 by Wojtek Pych, CAMK PAN           */
/*  pych@camk.edu.pl                                      */
/*                                                        */
/*  Written for GNU project C and C++ Compiler            */
/*                                                        */
/* Program is finding stars in the image and outputs      */
/* positions, magnitudes and background levels.           */
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
#include <stdlib.h>
#include <string.h>
#include <strings.h>

#include "defs.h"

#include "errmess.h"
#include "pfitsin.h"
#include "bkg.h"
#include "covar.h"
#include "indexx.h"
#include "centroid.h"
#include "aperphot.h"

void usage();
PARAMS read_params(char *, char *);

/*--------------------------------------------------------*/
void usage()
{
  printf("\n\tUSAGE: sfind  parameter_file instrument_file ");
  printf("input_image output_file\n\n");
  exit(EXIT_FAILURE);

  return;
}
/*--------------------------------------------------------*/
PARAMS read_params(char *parfname, char *instrname)
{
        char    line[201],
                key[200],
                val[200],
                algorithm[20];
        int     i;
        FILE    *inf;
        PARAMS  lpar;

  if (!(inf = fopen(parfname, "r"))) errmess(parfname);

  for (i=0; !feof(inf); i++)
  {
    if (fgets(line, 200, inf) == NULL)
    {
      printf("\n\tERROR! reading %s\n", parfname);
      exit(EXIT_FAILURE);
    }
    if (feof(inf)) break;

    sscanf(line, "%s = %s", key, val);

    if (!strcmp(key, "END")) break;
    else if (!strcmp(key, "MOHW"))        sscanf(val, "%d", &lpar.mohw);
    else if (!strcmp(key, "SMHW"))        sscanf(val, "%d", &lpar.smhw);
    else if (!strcmp(key, "N_SIG"))       sscanf(val, "%d", &lpar.n_sig);
    else if (!strcmp(key, "SIG_X"))       sscanf(val, "%lf", &lpar.sig_x);
    else if (!strcmp(key, "SIG_Y"))       sscanf(val, "%lf", &lpar.sig_y);
    else if (!strcmp(key, "C_MIN"))       sscanf(val, "%lf", &lpar.c_min);
    else if (!strcmp(key, "APRAD"))       sscanf(val, "%lf", &lpar.aprad);
    else if (!strcmp(key, "ANRAD1"))      sscanf(val, "%lf", &lpar.anrad1);
    else if (!strcmp(key, "ANRAD2"))      sscanf(val, "%lf", &lpar.anrad2);
    else if (!strcmp(key, "SIG_THRESH"))  sscanf(val, "%lf", &lpar.sig_thresh);
    else if (!strcmp(key, "ABS_THRESH"))  sscanf(val, "%lf", &lpar.abs_thresh);
    else if (!strcmp(key, "BKG_ALG"))     sscanf(val, "%s", algorithm);
    else if (!strcmp(key, "VERBOSE"))     sscanf(val, "%hd", &lpar.verbose);
  }

  fclose(inf);

  if (!strcasecmp(algorithm, "mode"))         lpar.bkg_algorithm=1;
  else if (!strcasecmp(algorithm, "median"))  lpar.bkg_algorithm=2;
  else if (!strcasecmp(algorithm, "mean"))    lpar.bkg_algorithm=3;
  else                                        lpar.bkg_algorithm=0;

  if (lpar.verbose > 1)
  {
    printf("Reading %s\n", parfname);
    printf("----------\n");
    printf("MOHW=       %d\n", lpar.mohw);
    printf("SMHW=       %d\n", lpar.smhw);
    printf("N_SIG=      %d\n", lpar.n_sig);
    printf("SIG_X=      %g\n", lpar.sig_x);
    printf("SIG_Y=      %g\n", lpar.sig_y);
    printf("C_MIN=      %g\n", lpar.c_min);
    printf("APRAD=      %g\n", lpar.aprad);
    printf("ANRAD1=     %g\n", lpar.anrad1);
    printf("ANRAD2=     %g\n", lpar.anrad2);
    printf("SIG_THRESH= %g\n", lpar.sig_thresh);
    printf("ABS_THRESH= %g\n", lpar.abs_thresh);
    printf("BKG_ALG=    %d\n", lpar.bkg_algorithm);
    printf("VERBOSE=    %d\n", lpar.verbose);
    printf("----------\n");
  }

  if (!(inf = fopen(instrname, "r"))) errmess(instrname);

  for (i=0; !feof(inf); i++)
  {
    if (fgets(line, 200, inf) == NULL)
    {
      printf("\n\tERROR! reading %s\n", instrname);
      exit(EXIT_FAILURE);
    }
    if (feof(inf)) break;

    sscanf(line, "%s = %s", key, val);

    if (!strcasecmp(key, "END")) break;
    else if (!strcasecmp(key, "SAT_LEVEL"))  sscanf(val, "%lf", &lpar.sat_level);
    else if (!strcasecmp(key, "MIN_LEVEL"))  sscanf(val, "%lf", &lpar.min_level);
    else if (!strcasecmp(key, "GAIN"))       continue;
    else if (!strcasecmp(key, "NX"))         continue;
    else if (!strcasecmp(key, "NY"))         continue;
    else if (!strcasecmp(key, "MARG"))       continue;
    else if (!strcasecmp(key, "FITS"))       continue;
    else printf("WARNING: Unknown parameter %s\n", key);
  }

  fclose(inf);

  if (lpar.verbose > 1)
  {
    printf("Reading %s\n", instrname);
    printf("----------\n");
    printf("SAT_LEVEL     = %g\n", lpar.sat_level);
    printf("MIN_LEVEL     = %g\n", lpar.min_level);
    printf("----------\n");
  }

  return(lpar);
}
/*--------------------------------------------------------*/
int max_test(int nx, int ny, double **corrim, int i, int j, int srad)
{
        int   k, l;
        double cpix;

  cpix=corrim[j][i];

  for (l=-srad; l<=srad; l++)
    if ((i+l >= 0) && (i+l < ny))
      for (k=-srad; k<=srad; k++)
        if ((j+k >= 0) && (j+k < nx))
          if (cpix < corrim[j+k][i+l]) return(0);

  return(1);
}
/*--------------------------------------------------------*/
double bckg_sigma(int nx, int ny, double **data, PARAMS par, double *sigma)
{
        int   lx, ly,         // subframe size
              i, j,
              m, n;
        double bkg[9],
              mbkg,
              stddev,
              min_stddev,
              **sdata;

  min_stddev=par.sat_level;
  lx=nx/3;
  ly=ny/3;

  if (!(sdata=(double **)calloc((size_t)ly, sizeof(double *))))
    errmess("bckg_sigma() calloc(sdata)");
  for (i=0; i<ly; i++)
    if (!(sdata[i]=(double *)calloc((size_t)lx, sizeof(double))))
      errmess("bckg_sigma() calloc(sdata[i])");

  for (i=0; i<3; i++)
  {
    for (j=0; j<3; j++)
    {
      for (m=0; m<ly; m++)
      {
        for (n=0; n<lx; n++)
        {
          sdata[m][n]=data[i*ly+m][j*lx+n];
        }
      }
      bkg[i*3+j]=bkg_level(lx, ly, sdata, par, &stddev);
      if (stddev < min_stddev) min_stddev=stddev;
    }
  }

  *sigma=min_stddev;
  mbkg=0.0;
  for (i=0; i<9; i++) mbkg+=bkg[i];
  mbkg/=9.0;

  return(mbkg);
}
/*--------------------------------------------------------*/
int main(int argc, char *argv[])
{
        char    **header,
                *imname,
                *parfname,
                *instrname,
                *outfname;
        int     nx,
                ny,
                i, j,
                marg,         // margin from border for correct position
                sat_test,
                nstar,
                headlen;
        double   **im,
                **corrim,
                x,
                y,
                bglevel,
                stdev,
                phot,
                thresh;
        FILE    *outf;
        PARAMS  par;

/* IO stuff */
  if (argc != 5) usage();

  parfname = argv[1];
  instrname= argv[2];
  imname   = argv[3];
  outfname = argv[4];

  par=read_params(parfname, instrname);

  im=read_FITS_2Dfile(imname, 's', &headlen, &header, &nx, &ny);
  if (par.verbose > 1) printf("%s: %d x %d [pixels]\n", imname, nx, ny);

  for (i=0; i<headlen; i++) free(header[i]);
  free(header);

/* find background and std_dev for threshold */
  thresh=bckg_sigma(nx, ny, im, par, &stdev);
  if (par.verbose > 1) printf("background= %g\n", thresh);

  if (par.abs_thresh < 0.0) thresh+=par.sig_thresh*stdev;
  else                      thresh+=par.abs_thresh;
  if (par.verbose > 1) printf(" detection threshold= %g\n", thresh);

/* calculate covariance matrix by convolving with gaussian model */
  corrim=covar(nx, ny, im, thresh, par);

  nstar=0;

  if (!(outf=fopen(outfname, "w"))) errmess(outfname);

/* look for local maxima in covariance matrix */
  marg=(int)(par.aprad/2.0);
  if (par.verbose > 1) printf(" margin= %d\n", marg);

  for (j=marg; j<ny-marg; j++)
  {
    for (i=marg; i<nx-marg; i++)
    {
      if (im[j][i] > thresh)
      {
        if (max_test(nx, ny, corrim, i, j, par.smhw))
        {
          centroid(nx, ny, im, i, j, &x, &y);

          if (corrim[j][i] > par.c_min)
          {
            bglevel=bkg(nx, ny, im, i, j, par.anrad1, par.anrad2, par);
            phot=aperphot(nx, ny, im, i, j, bglevel, &sat_test, par);
            fprintf(outf, "%10.4f %10.4f %9.3f %9.3f  %d\n",
                             x, y, phot, bglevel, sat_test);
            nstar++;
          }
        }
      }
    }
  }

  fclose(outf);

  if (par.verbose) printf("%d stars found\n", nstar);

  for (i=0; i<ny; i++) free(im[i]);
  free(im);
  free(corrim);

  return(0);
}
/*** END ***/
