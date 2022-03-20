/*========================================================*/
/*                                                        */
/*  bkg.c               version 1.2.3   2005.01.24        */
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
/* Given image im(nx, ny) and centroid (cx, cy)           */
/* estimates background taking data from annulus (r1, r2) */
/* Sigma clipping with n_sig until nothing changes        */
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

#include <stdlib.h>

#include "defs.h"
#include "funcs.h"
#include "errmess.h"

double bkg(double *im, int cx, int cy, PAR_STRUCT *par)
{
        int   i, j, npix, ir, *idx, delta, nx, ny;
        double r1, r2, rr, rr1, rr2, bglev, absdev, mean, *tab, n_sig;

  n_sig = par->nsig_bkg;

  r1 = par->anrad1;
  r2 = par->anrad2;

  nx = par->nx0;
  ny = par->ny0;

  npix = 3.1416*((r2+1.0)*(r2+1.0) - r1*r1);

  if (!(tab  = (double *) malloc(npix*sizeof(double)))) errmess("malloc(tab)");
  if (!(idx  =   (int *) malloc(npix*sizeof(int)))) errmess("malloc(tab)");

  ir   = r2 + 1;
  rr1  = r1*r1;
  rr2  = r2*r2;

/* load pixels within annulus (r1, r2) into tab */
  npix = 0;

  for (j=-ir; j<=ir; j++)
  {
    for (i=-ir; i<=ir; i++)
    {
      rr = i*i + j*j;

      if (rr1<rr && rr<rr2) tab[npix++] = im[cx+i+nx*(cy+j)];
    }
  }

/* take median as first background level                             */
/* use little trick to fool stupid 1,n indexing of numerical recipes */

  indexx(npix, tab-1, idx-1);

  bglev = tab[idx[npix/2]-1];

/* reject pixels deviating by more than n_sig*sigma from current bglev */
/* subsequent bglev will be means of the remaining pixels              */
/* quit when no more rejections                                        */

  delta = 1;

  while (delta > 0)
  {
    absdev=0.0;

    for (i=0; i<npix; i++)
      absdev += (bglev - tab[i])*(bglev - tab[i]);

    absdev *= n_sig*n_sig/(double)npix;

    j=0;
    mean=0.0;

    for (i=0; i<npix; i++)
    {
      if ((bglev - tab[i])*(bglev - tab[i]) < absdev)
      {
        tab[j++] = tab[i];
        mean += tab[i];
      }
    }

    if (j > 1)
    {
      delta = npix - j;
      npix  = j - 1;
      mean /= npix;

      bglev = mean;
    }
    else  return(0.0);
  }

  free(tab);
  free(idx);

  return(bglev);
}
