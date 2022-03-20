/*========================================================*/
/*                                                        */
/*  xymatch.c           version 1.4.5   2006.10.23        */
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
/* Program matches two (x,y) cordinate lists              */
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
#include <string.h>

#include "defs.h"
#include "errmess.h"

#include "triangles.h"
#include "bright_end.h"
#include "xy_lin.h"
#include "refine.h"

void usage();
void read_params(char *, PARAMS *);
int  read_objects(char *, double **, double **, double **, short **);

/*--------------------------------------------------------*/
void usage()
{
  printf("\n\t USAGE: xymatch  parameter_file ");
  printf("reference_list matched_list output_file corr_file_name\n\n");
  exit(1);
}
/*--------------------------------------------------------*/
void read_params(char *iname, PARAMS *par)
{
        char  line[201],
              key[200],
              val[200];
        int   i;
        FILE  *inf;

  if (!(inf=fopen(iname, "r"))) errmess(iname);
  for (i=0; !feof(inf); i++)
  {
    fgets(line, 200, inf);
    if (feof(inf)) break;

    sscanf(line, "%s = %s", key, val);

    if (!strcasecmp(key, "END")) break;
    else if (!strcasecmp(key, "NSUB"))    sscanf(val, "%d", &par->nsub);
    else if (!strcasecmp(key, "LLIM"))    sscanf(val, "%lf", &par->llim);
    else if (!strcasecmp(key, "RLIM"))    sscanf(val, "%lf", &par->rlim);
    else if (!strcasecmp(key, "FVNO"))    sscanf(val, "%lf", &par->fvno);
    else if (!strcasecmp(key, "LTOL"))    sscanf(val, "%lf", &par->ltol);
    else if (!strcasecmp(key, "RTOL"))    sscanf(val, "%lf", &par->rtol);
    else if (!strcasecmp(key, "CTOL"))    sscanf(val, "%lf", &par->ctol);
    else if (!strcasecmp(key, "PTOL"))    sscanf(val, "%lf", &par->ptol);
    else if (!strcasecmp(key, "VERBOSE")) sscanf(val, "%hd", &par->verbose);
    else printf("WARNING: Unknown parameter %s\n", key);
  }

  fclose(inf);

  if (par->verbose > 1)
  {
    printf("Reading %s:\n", iname);
    printf("-----------\n");
    printf("NSUB = %d\n", par->nsub);
    printf("LLIM = %g\n", par->llim);
    printf("RLIM = %g\n", par->rlim);
    printf("FVNO = %g\n", par->fvno);
    printf("LTOL = %g\n", par->ltol);
    printf("RTOL = %g\n", par->rtol);
    printf("CTOL = %g\n", par->ctol);
    printf("PTOL = %g\n", par->ptol);
    printf("VERBOSE = %hd\n", par->verbose);
    printf("-----------\n");
  }

  return;
}
/*--------------------------------------------------------*/
int read_objects(char *iname, double **x, double **y, double **mag, short **sat)
{
        char  buf[256];
        short *ls;
        int   i;
        double *lx, *ly, *lm;
        FILE  *inf;

  if (!(lx=(double *)calloc(1, sizeof(double)))) errmess("calloc(lx)");
  if (!(ly=(double *)calloc(1, sizeof(double)))) errmess("calloc(ly)");
  if (!(lm=(double *)calloc(1, sizeof(double)))) errmess("calloc(lm)");
  if (!(ls=(short *)calloc(1, sizeof(short)))) errmess("calloc(ls)");

  if (!(inf=fopen(iname, "r"))) errmess(iname);

  for (i=0; !feof(inf); i++)
  {
    fgets(buf, 256, inf);
    if (feof(inf)) break;
    
    if (!(lx=(double *)realloc(lx, (i+1)*sizeof(double)))) errmess("realloc(lx)");
    if (!(ly=(double *)realloc(ly, (i+1)*sizeof(double)))) errmess("realloc(ly)");
    if (!(lm=(double *)realloc(lm, (i+1)*sizeof(double)))) errmess("realloc(lm)");
    if (!(ls=(short *)realloc(ls, (i+1)*sizeof(short)))) errmess("realloc(ls)");

    sscanf(buf, "%lf %lf %lf %*s %hd", &lx[i], &ly[i], &lm[i], &ls[i]);

    if (lm[i] < -999.99) i--; // this means wrong photometry from sfind
  }

  fclose(inf);

  *x=lx;
  *y=ly;
  *mag=lm;
  *sat=ls;

  return(i);
}
/*--------------------------------------------------------*/
int reject_saturated(int nobj, double **x, double **y, short *sat)
{
        int   i,
              nn;
        double *nx, *ny,
              *lx, *ly;

  lx=*x;
  ly=*y;

  if (!(nx=(double *)calloc(1, sizeof(double))))
    errmess("reject_saturated: calloc(nx)");
  if (!(ny=(double *)calloc(1, sizeof(double))))
    errmess("reject_saturated: calloc(ny)");

  nn=0;
  for (i=0; i<nobj; i++)
  {
    if (!sat[i])
    {
      if (!(nx=(double *)realloc(nx, (nn+1)*sizeof(double))))
        errmess("reject_saturated: realloc(nx)");
      if (!(ny=(double *)realloc(ny, (nn+1)*sizeof(double))))
        errmess("reject_saturated: realloc(ny)");

      nx[nn]=lx[i];
      ny[nn]=ly[i];
      nn++;
    }
  }

  free(lx);
  free(ly);

  *x=nx;
  *y=ny;

  return(nn);
}
/*--------------------------------------------------------*/
int main(int argc, char *argv[])
{
        char    *inp1name, *inp2name, *parfname, *outfname,
                *corrname;
        short   *sat1, *sat2;
        int     i,
                nobj1, nobj2,
                nsub1, nsub2,
                nmatch,
                maxobj,
                *index;
        double   *x1, *y1, *x2, *y2,
                *mag1, *mag2,
                *xs1, *ys1, *xs2, *ys2,
                *xm1, *ym1, *xm2, *ym2,
                coeffx[3], coeffy[3];
        FILE    *outf;
        PARAMS  par;

/* IO stuff */
  if (argc != 6) usage();

  parfname = argv[1];
  inp1name = argv[2];
  inp2name = argv[3];
  outfname = argv[4];
  corrname = argv[5];

  read_params(parfname, &par);

  nobj1=read_objects(inp1name, &x1, &y1, &mag1, &sat1);
  if (par.verbose) printf("%d objects read from %s\n", nobj1, inp1name);
  if (par.nsub > nobj1) par.nsub = nobj1;
  nsub1=bright_end(nobj1, x1, y1, mag1, par.nsub, &xs1, &ys1, par.verbose);

  nobj2=read_objects(inp2name, &x2, &y2, &mag2, &sat2);
  if (par.verbose) printf("%d objects read from %s\n", nobj2, inp2name);
  if (par.nsub > nobj2) par.nsub = nobj2;
  nsub2=bright_end(nobj2, x2, y2, mag2, par.nsub, &xs2, &ys2, par.verbose);

  par.nsub=(nsub1 < nsub2 ? nsub1 : nsub2);
  if (par.verbose > 1) printf("par.nsub= %d\n", par.nsub);

  free(mag1);
  free(mag2);

  if (par.verbose > 2)
  {
    printf("Coordinates of the brighest objects:\n");
    printf("  X1    Y1      X2    Y2\n");
    printf("--------------------------\n");
    for (i=0; i<par.nsub; i++)
      printf("%8.2f %8.2f    %8.2f %8.2f\n", xs1[i], ys1[i], xs2[i], ys2[i]);
    printf("--------------------------\n\n");
  }

/* match nsub brightest stars for approximate transformation */
  maxobj=(nobj1 > nobj2 ? nobj1 : nobj2);
  if (par.verbose > 2) printf("maxobj= %d\n", maxobj);

  if (!(index=(int *)calloc(maxobj, sizeof(int)))) errmess("calloc(index)");

  triangles(xs1, ys1, xs2, ys2, par.nsub, par.nsub, index, par);

  if (!(xm1=(double *)malloc(sizeof(double)))) errmess("malloc(xm1)");
  if (!(ym1=(double *)malloc(sizeof(double)))) errmess("malloc(ym1)");
  if (!(xm2=(double *)malloc(sizeof(double)))) errmess("malloc(xm2)");
  if (!(ym2=(double *)malloc(sizeof(double)))) errmess("malloc(ym2)");

  nmatch=0;
  for (i=0; i<par.nsub; i++)
  {
    if (index[i] != -1)
    {
      if (!(xm1=(double *)realloc(xm1, (nmatch+1)*sizeof(double))))
        errmess("realloc(xm1)");
      if (!(ym1=(double *)realloc(ym1, (nmatch+1)*sizeof(double))))
        errmess("realloc(ym1)");
      if (!(xm2=(double *)realloc(xm2, (nmatch+1)*sizeof(double))))
        errmess("realloc(xm2)");
      if (!(ym2=(double *)realloc(ym2, (nmatch+1)*sizeof(double))))
        errmess("realloc(ym2)");

      xm1[nmatch]=xs1[i];
      ym1[nmatch]=ys1[i];
      xm2[nmatch]=xs2[index[i]];
      ym2[nmatch]=ys2[index[i]];
      nmatch++;
    }
  }

  free(xs1);
  free(ys1);
  free(xs2);
  free(ys2);

  if (nmatch < 2)
  {
    printf("ERROR: nmatch < 2\n");
    exit(2);
  }
  if (par.verbose) printf("%d objects matched by triangles()\n", nmatch);

/* linear fit to nmatch stars indentified by triangles */
  xy_lin(xm1, ym1, xm2, ym2, nmatch, coeffx, coeffy);

  free(xm1);
  free(ym1);
  free(xm2);
  free(ym2);

  if (par.verbose > 1)
  {
    printf("Linear transformation:\n");
    printf("----------------------\n");
    for (i=0; i<3; i++)
      printf("coeffx[%d]= %12g   coeffy[%d]= %12g\n",
        i, coeffx[i], i, coeffy[i]);
    printf("----------------------\n");
  }

  nobj1=reject_saturated(nobj1, &x1, &y1, sat1);
  nobj2=reject_saturated(nobj2, &x2, &y2, sat2);
  if (par.verbose > 1)
  {
    printf("%d objects from %s left after reject_saturated()\n",
      nobj1, inp1name);
    printf("%d objects from %s left after reject_saturated()\n",
      nobj2, inp2name);
  }

  free(sat1);
  free(sat2);

/* using linear fit transform one list and look for close neighbors */
  for (i=0; i<nobj1; i++)
    maxobj=refine(coeffx, coeffy, x1, y1, x2, y2, nobj1, nobj2, index,
      par.ptol);
  if (par.verbose)
  {
    printf("%d objects left in the template list after refine()\n", maxobj);
    printf("Writing matched list to %s\n\n", outfname);
  }

  if (!(outf=fopen(outfname, "w"))) errmess("outfname");
  for (i=0; i<nobj1; i++)
    if (index[i] != -1)
      fprintf(outf, "%9.3f %10.3f %10.3f %10.3f\n",
                    x1[i], y1[i], x2[index[i]], y2[index[i]]);
  fclose(outf);

  free(index);
  free(x1); free(y1);
  free(x2); free(y2);

  if (!(outf=fopen(corrname, "w"))) errmess(corrname);

  fprintf(outf, "%d  %d  %s",
          (int)(coeffx[2]+0.5), (int)(coeffy[2]+0.5), inp2name);

  fclose(outf);

  return(0);
}
/*** END ***/
