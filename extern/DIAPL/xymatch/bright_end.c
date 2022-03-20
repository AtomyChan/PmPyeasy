/*========================================================*/
/*                                                        */
/*  bright_end.c        version 1.4.1   2006.03.31        */
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
/* Given a list of nd objects:                            */
/* (position, magnitude) = (x, y, mag)                    */
/* return a similar list for nobj brightest objects.      */
/* Objects are spread in 9 subrasters.                    */
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

#include "indexx.h"
#include "errmess.h"

/*--------------------------------------------------------*/
int subselect(int nd, float *x, float *y, float *mag,
              float x1, float x2, float y1, float y2,
              float **tx, float **ty, float **tmag)
{
        int   i,
              ns;
        float *lx, *ly, *lmag;

  if (!(lx=(float *)calloc(nd, sizeof(float))))
    errmess("subselect(): calloc(lx)");
  if (!(ly=(float *)calloc(nd, sizeof(float))))
    errmess("subselect(): calloc(ly)");
  if (!(lmag=(float *)calloc(nd, sizeof(float))))
    errmess("subselect(): calloc(lmag)");

  ns=0;
  for (i=0; i<nd; i++)
  {
    if ((x[i] > x1) && (x[i] < x2) && (y[i] > y1) && (y[i] < y2))
    {
      lx[ns]=x[i];
      ly[ns]=y[i];
      lmag[ns]=mag[i];
      ns++;
    }
  }

  if (!(lx=(float *)realloc(lx, ns*sizeof(float))))
    errmess("subselect(): realloc(lx)");
  if (!(ly=(float *)realloc(ly, ns*sizeof(float))))
    errmess("subselect(): realloc(ly)");
  if (!(lmag=(float *)realloc(lmag, ns*sizeof(float))))
    errmess("subselect(): realloc(lmag)");

  *tx=lx;
  *ty=ly;
  *tmag=lmag;

  return(ns);
}
/*--------------------------------------------------------*/
int bright_end(int nobj, float *x, float *y, float *mag,
                int nbright, float **xs, float **ys, short verbose)
{
        int   i, j, k,
              nmax,
              m,
              subnbright,
              nsel,
              *idx;
        float *lxs, *lys,
              *tx, *ty, *tmag,
              xmin, xmax, ymin, ymax,
              x1, x2, y1, y2,
              dx, dy;

  if (nbright < 9) nbright=9;
  subnbright=nbright/9;
  nbright=subnbright*9;

  if (verbose > 1)
    printf("bright_end():\n\t subnbright= %d  nbright= %d  nobj= %d\n",
      subnbright, nbright, nobj);

  if (nbright > nobj)
  {
    printf("ERROR! bright_end(): not enough objects in the list\n");
    exit(3);
  }

/* check coordinates range */
  xmin = xmax = x[0];
  ymin = ymax = y[0];
  for (i=1; i<nobj; i++)
  {
    if (x[i] < xmin) xmin=x[i];
    if (x[i] > xmax) xmax=x[i];
    if (y[i] < ymin) ymin=y[i];
    if (y[i] > ymax) ymax=y[i];
  }

  dx=(xmax-xmin)/3.0;
  dy=(ymax-ymin)/3.0;
  if (verbose > 1)
  {
    printf("\t coordinates range: [%.2f : %.2f, %.2f : %.2f]\n",
      xmin, xmax, ymin, ymax);
    printf("\t dx= %g  dy= %g\n", dx, dy);
  }

  if (!(lxs=(float *)calloc(nbright, sizeof(float)))) errmess("calloc(lxs)");
  if (!(lys=(float *)calloc(nbright, sizeof(float)))) errmess("calloc(lys)");

  m=0;
  for (i=0; i<3; i++)
  {
    for (j=0; j<3; j++)
    {
      x1=xmin+i*dx; x2=xmin+(i+1)*dx;
      y1=ymin+j*dy; y2=ymin+(j+1)*dy;
      nsel=subselect(nobj, x, y, mag, x1, x2, y1, y2, &tx, &ty, &tmag);
      if (verbose > 2)
        printf("[%.2f : %.2f, %.2f : %.2f]  nsel= %d\n", x1, x2, y1, y2, nsel);
      if (!nsel) continue;
      if (nsel == 1)
      {
        lxs[m]=tx[0];
        lys[m]=ty[0];

        m++;

        free(tx);
        free(ty);
        free(tmag);

        continue;
      }

      if (!(idx=(int *)malloc(nsel*sizeof(int)))) errmess("malloc(idx)");

/* use little trick to fool stupid 1,n indexing of numerical recipes */
      indexx(nsel, tmag-1, idx-1);

      nmax=(nsel < subnbright ? nsel : subnbright);
      if (verbose > 2) printf("\t nmax= %d\n", nmax);

      for (k=0; k<nmax; k++)
      {
        lxs[m+k]=tx[idx[k]-1];
        lys[m+k]=ty[idx[k]-1];
      }
      m+=nmax;

      free(tx);
      free(ty);
      free(tmag);
      free(idx);
    }
  }

  *xs=lxs;
  *ys=lys;

  return(m);
}
