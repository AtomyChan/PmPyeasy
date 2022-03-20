/*========================================================*/
/*                                                        */
/*  xygrid.c            version 1.1.6   2005.04.07        */
/*                                                        */
/*  Original source code:                                 */
/*  Copyright (C) 2005 by Przemek Wozniak                 */
/*  wozniak@lanl.gov                                      */
/*                                                        */
/*  Modifications:                                        */
/*  Copyright (C) 2005 by Wojtek Pych, CAMK PAN           */
/*  pych@camk.edu.pl                                      */
/*                                                        */
/*  Written for GNU project C and C++ Compiler            */
/*                                                        */
/* Given matched coordinate lists of objects identified   */
/* in two frames calculates the coordinate transformation */
/* of the specified order.                                */
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

#include <math.h>
#include <stdio.h>
#include <string.h>

#include "fitn.h"
#include "poly.h"
#include "sigma.h"
#include "errmess.h"

typedef struct
{
  short verbose;
  int   ndeg,
        maxniter;
  float sigmaf;
} PARAMS;

/*--------------------------------------------------------*/
void readpar(char *parfname, PARAMS *par)
{
        char  line[201],
              key[200],
              val[200];
        FILE  *inpf;

  if (!(inpf = fopen(parfname, "r"))) errmess(parfname);

  while (!feof(inpf))
  {
    fgets(line, 200, inpf);
    if (feof(inpf)) break;

    sscanf(line, "%s = %s", key, val);

    if (!strcasecmp(key, "END")) break;
    else if (!strcasecmp(key, "NDEG"))      sscanf(val, "%d", &par->ndeg);
    else if (!strcasecmp(key, "MAX_NITER")) sscanf(val, "%d", &par->maxniter);
    else if (!strcasecmp(key, "SIGMA_F"))   sscanf(val, "%f", &par->sigmaf);
    else if (!strcasecmp(key, "VERBOSE"))   sscanf(val, "%hd", &par->verbose);
    else printf("WARNING: Unknown parameter %s\n", key);
  }

  fclose(inpf);

  if (par->verbose > 1)
  {
    printf("Reading %s\n", parfname);
    printf("----------\n");
    printf("NDEG      = %d\n", par->ndeg);
    printf("MAX_NITER = %d\n", par->maxniter);
    printf("SIGMA_F   = %g\n", par->sigmaf);
    printf("VERBOSE   = %d\n", par->verbose);
    printf("----------\n");
  }

  return;
}
/*--------------------------------------------------------*/
int readcoor(char *inpfname, float **x, float **y, float **zx, float **zy)
{
        char  buf[256];
        int   i;
        float *lx, *ly, *lzx, *lzy;
        FILE  *inpf;

  if (!(lx=(float *)calloc(1, sizeof(float))))
    errmess("readcoor: calloc(lx)");
  if (!(ly=(float *)calloc(1, sizeof(float))))
    errmess("readcoor: calloc(ly)");
  if (!(lzx=(float *)calloc(1, sizeof(float))))
    errmess("readcoor: calloc(lzx)");
  if (!(lzy=(float *)calloc(1, sizeof(float))))
    errmess("readcoor: calloc(lzy)");

  if (!(inpf=fopen(inpfname, "r"))) errmess(inpfname);

  for (i=0; !feof(inpf); i++)
  {
    fgets(buf, 256, inpf);
    if (feof(inpf)) break;

    if (!(lx=(float *)realloc(lx, (i+1)*sizeof(float))))
      errmess("readcoor: realloc(lx)");
    if (!(ly=(float *)realloc(ly, (i+1)*sizeof(float))))
      errmess("readcoor: realloc(ly)");
    if (!(lzx=(float *)realloc(lzx, (i+1)*sizeof(float))))
      errmess("readcoor: realloc(lzx)");
    if (!(lzy=(float *)realloc(lzy, (i+1)*sizeof(float))))
      errmess("readcoor: realloc(lzy)");

    sscanf(buf, "%f %f %f %f", &lx[i], &ly[i], &lzx[i], &lzy[i]);
  }

  fclose(inpf);

  *x=lx;
  *y=ly;
  *zx=lzx;
  *zy=lzy;

  return(i);
}
/*--------------------------------------------------------*/
int main(int argc, char *argv[])
{	
        char    *parfname, *inpfname, *outfname;
        int     i, nobj, nterm, ik, ik0, delta, niter;
        float   *x, *nx, *nzx, *zx,
                *y, *ny, *nzy, *zy,
                dx, dy, rr, thresh, sx, sy;
        double  *coeffx, *coeffy;
        FILE    *outf;
        PARAMS  par;

/* IO stuff */

  if (argc != 4)
  {
    printf("\n\tUSAGE: xygrid  parameter_file  input_list  coeff_file\n\n");
    exit(1);
  }

  parfname = argv[1];
  inpfname = argv[2];
  outfname = argv[3];

  readpar(parfname, &par);

  nterm = (par.ndeg+1)*(par.ndeg+2)/2;

  nobj=readcoor(inpfname, &x, &y, &zx, &zy);

  if (!(coeffx=(double *)calloc(nterm, sizeof(double))))
    errmess("calloc(coeffx)");
  if (!(coeffy=(double *)calloc(nterm, sizeof(double))))
    errmess("calloc(coeffy)");

/* make initial fit */
  fit(x, y, zx, par.ndeg, nobj, coeffx);
  fit(x, y, zy, par.ndeg, nobj, coeffy);

  sx = sigma(x, y, zx, par.ndeg, nobj, coeffx);
  sy = sigma(x, y, zy, par.ndeg, nobj, coeffy);

  thresh = par.sigmaf*par.sigmaf*(sx*sx + sy*sy);

  ik = nobj;
  delta = 1;
  niter = 0;

  if (!(nx=(float *)calloc(nobj, sizeof(float))))
    errmess("readcoor: calloc(nx)");
  if (!(ny=(float *)calloc(nobj, sizeof(float))))
    errmess("readcoor: calloc(ny)");
  if (!(nzx=(float *)calloc(nobj, sizeof(float))))
    errmess("readcoor: calloc(nzx)");
  if (!(nzy=(float *)calloc(nobj, sizeof(float))))
    errmess("readcoor: calloc(nzy)");

/* sigma clipping of the fit until MAX_NITER reached or nothing changes */

  while ((delta > 0) && (niter < par.maxniter))
  {
    ik0 = ik;
    ik  = 0;

    for (i=0; i<nobj; i++)
    {
      dx  = poly(x[i], y[i], par.ndeg, coeffx) - zx[i];
      dy  = poly(x[i], y[i], par.ndeg, coeffy) - zy[i];

      rr  = dx*dx + dy*dy;

      nx[ik]  = x[i];
      ny[ik]  = y[i];
      nzx[ik] = zx[i];
      nzy[ik] = zy[i];

      if (rr < thresh) ++ik;
    }

    delta = ik0 - ik;

    fit(nx, ny, nzx, par.ndeg, ik, coeffx);
    fit(nx, ny, nzy, par.ndeg, ik, coeffy);

    sx = sigma(nx, ny, nzx, par.ndeg, ik, coeffx);
    sy = sigma(nx, ny, nzy, par.ndeg, ik, coeffy);

    niter++;
  }

  free(x);
  free(y);
  free(zx);
  free(zy);
  free(nx);
  free(ny);
  free(nzx);
  free(nzy);

/* print results and store coefficients in binary file */

  if (par.verbose)
  {
    printf("\n");
    for (i=0; i<nterm; i++)
    {
      printf("coeffx[%d] = %9.6f   ", i, coeffx[i]);
      printf("coeffy[%d] = %9.6f \n", i, coeffy[i]);
    }
  }

  printf("%s:  sigmax= %.4f   sigmay= %.4f   ndata= %d   nleft= %d\n",
          outfname, sx, sy, nobj, ik);

  if (!(outf=fopen(outfname, "w"))) errmess(outfname);

  fwrite(&nterm, sizeof(int), 1, outf);
  fwrite(coeffx, sizeof(double), nterm, outf);
  fwrite(coeffy, sizeof(double), nterm, outf);

  fclose(outf);

  free(coeffx);
  free(coeffy);

  return(0);
}
/*** END ***/
