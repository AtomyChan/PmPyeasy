/*========================================================*/
/*                                                        */
/*  covar_sig.c         version 1.2.3   2005.01.24        */
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
/* Convolves image im(nx, ny)                             */
/* with the PSF model (2*hm+1)x(2*hm+1) from psfim,       */
/* corrim on output equals to covariance (signed)         */
/* of the image and model or 0.0                          */
/* for pixels below thresh abs value.                     */
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

void covar_sig(double *difim, double *corrim, int nx, int ny, double thresh,
               int hm, double *psfim, int psfhw)
{
        int   i, j, k, l, m, psfn;
        double norm, mean_im, mean_model, *model, a0, a, b, c;

  m     = 2*hm+1;
  norm  = m*m;
  psfn  = 2*psfhw+1;

  if (!(model = (double *)malloc((int)norm*sizeof(double))))
    errmess("malloc(model)");

/* prepare lowered gaussian filter */

  mean_model = 0.0;

  for (j=-hm; j<=hm; j++)
  {
    for (i=-hm; i<=hm; i++)
    {
      model[hm+i+m*(hm+j)] = psfim[psfhw+i+psfn*(psfhw+j)];
      mean_model += model[hm+i+m*(hm+j)];
    }
  }

  mean_model /= norm;

  a = 0.0;
  for (j=-hm; j<=hm; j++)
  {
    for(i=-hm; i<=hm; i++)
    {
      a0 = model[hm+i+m*(hm+j)] - mean_model;
      a += a0*a0;
    }
  }

  for (j=0; j<ny; j++)
    for(i=0; i<nx; i++)
      corrim[i+nx*j] = 0.0;

/* convolve and normalize to get covariance */

  for (j=hm; j<ny-hm; j++)
  {
    for(i=hm; i<nx-hm; i++)
    {
      if (fabs(difim[i+nx*j]) > thresh)
      {
        b=c=0.0;
        mean_im=0.0;

        for (l=-hm; l<=hm; l++)
          for (k=-hm; k<=hm; k++)
            mean_im += difim[i+k+nx*(j+l)];

        mean_im /= norm;

        for (l=-hm; l<=hm; l++)
        {
          for (k=-hm; k<=hm; k++)
          {
            a0 = difim[i+k+nx*(j+l)] - mean_im;
            c+=a0*a0;
            b+=(difim[i+k+nx*(j+l)]-mean_im)*(model[hm+k+m*(hm+l)]-mean_model);
          }
        }

        if ((a != 0.0) && (c != 0.0)) corrim[i+nx*j]=b/sqrt(a*c);
        else                          corrim[i+nx*j]=0.0;
      }
      else  corrim[i+nx*j]=0.0;
    }
  }

  free(model);

  return;
}
