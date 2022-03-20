/*========================================================*/
/*                                                        */
/*  rm_bad_col.c        version 1.2.4   2005.02.16        */
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
/* Remove bad column. Every bad columns on image im       */
/* are replaced by gaussian noise.                        */
/* Bad column is a column, which mean luminosity          */
/* is greater than SAT_LEVEL or neighbor of bad column,   */
/* which mean luminosity is greater than 3*bg.            */
/* (C) Igor Soszynski & Karol Zebrun                      */
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
#include <stdlib.h>

#include "errmess.h"
#include "indexx.h"
#include "gasdev.h"

#include "defs.h"

double rm_bad_col(RINGPAR par, double *im1, int nx, int ny, int *col)
{
        int     i, j, n, npix, *idx, init;
        double   med, pix;
        double  mean, c;

  if (par.verbose > 1) printf("rm_bad_col()\n------------\n");

  npix = nx*ny;

  if (!(idx = (int *)malloc(npix*sizeof(int)))) errmess("malloc(idx)");

  indexx(npix, im1-1, idx-1);

  n = 0;
  while (((im1[idx[n]] == 0.0) || (im1[idx[n+1]] == 0.0)) && (n < npix)) n++;
  if (par.verbose > 2) printf("(im=0) => n/npix= %d/%d\n", n, npix);

  if (n != npix)  med = (double)im1[idx[(int)((npix-n)*0.3+n)]-1];
  else
  {
    if (par.verbose > 2) printf("Zero!\n");
    return(0.0);
  }

  c = sqrt(med/par.gain);
  if (par.verbose > 2) printf("med= %g   c= %g\n", med, c);

  for (i=0; i<nx; i++)
  {
    mean = 0.0;
    n = 0;

    for (j=0; j<ny; j++)
    {
      pix = im1[j*nx+i];
      if (pix > 0)
      {
        mean += pix;
        n++;
      }
    }

    if (n > 0)
    {
      mean/=n;
      col[i] = (int)mean;
    }
    else        col[i] = 0;

    if (par.verbose > 3)
      printf("column: %d  mean= %g  n= %d\n", i, mean, n);
  }

  if (par.growrad)
  {
    if (par.verbose > 1) printf("warning: growing radius\n");

    for (i=0; i<nx-1; i++)
      if (col[i+1] > par.sat)
        for (j=0; j<par.growrad; j++)
          if (i-j >= 0)
            col[i-j] = par.sat + 1;

    for (i=nx-1; i>0; i--)
      if (col[i-1] > par.sat)
        for (j=0; j<par.growrad; j++)
          if (i+j < nx)
            col[i+j] = par.sat + 1;
  }

  for (i=1; i<nx-1; i++)
    if (col[i] > par.nsigrm*med)
      if ((col[i-1] > par.sat) || (col[i+1] > par.sat))
        col[i] = par.sat + 1;

  n=0;
  for (i=0; i<nx; i++)
  {
    if (col[i] > par.sat)
    {
      if (par.verbose > 2) printf("column %d replaced\n", i);
      n++;
      for (j=0; j<ny; j++)
      {
        if (im1[j*nx+i] > 0.0)
        {
          init = i;
          im1[j*nx+i] = gasdev(&init)*c+med;
        }
      }
    }
  }

  if (par.verbose > 1) printf("%d columns replaced\n", n);

  return(med);
}
/*** END ***/
