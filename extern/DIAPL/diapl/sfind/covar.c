/*========================================================*/
/*                                                        */
/*  covar.c             version 1.2.5   2006.10.23        */
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
/* Convolves image im(nx,ny) with the gaussian model      */
/* (2*hm+1)x(2*hm+1) of sigmas (sig_x,sig_x).             */
/* corrim on output equals to covariance of the image     */
/* and model or 0.0 for pixels below thresh.              */
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

#include "defs.h"
#include "errmess.h"

double **covar(int nx, int ny, double **im, double thresh, PARAMS par)
{
        int   i, j, k, l, m,
              norm;
        double mean_im, a0, a, b, c,
              **gauss, mean_model,
              **corrim;

  m=2*par.mohw+1;

  if (!(gauss=(double **)calloc((size_t)m, sizeof(double))))
    errmess("calloc(gauss)");
  gauss+=par.mohw;
  for (j=-par.mohw; j<=par.mohw; j++)
  {
    if (!(gauss[j]=(double *)calloc((size_t)m, sizeof(double))))
      errmess("calloc(gauss[j])");
    gauss[j]+=par.mohw;
  }

/* prepare lowered gaussian filter */
  mean_model=0.0;

  for (j=-par.mohw; j<=par.mohw; j++)
  {
    for (i=-par.mohw; i<=par.mohw; i++)
    {
      a=(double)i/par.sig_x;
      b=(double)j/par.sig_y;
      gauss[j][i]=exp(-(a*a+b*b)/2.0);
      mean_model+=gauss[j][i];
    }
  }
  mean_model/=m*m;

  a=0.0;
  for (j=-par.mohw; j<=par.mohw; j++)
  {
    for(i=-par.mohw; i<=par.mohw; i++)
    {
      a0=gauss[j][i]-mean_model;
      a+=a0*a0;
    }
  }
  if (a <= 0.0)
  {
    printf("ERROR! covar: a <= 0.0 - this shoud not happen\n");
    exit(EXIT_FAILURE);
  }

  if (!(corrim=(double **)calloc((size_t)ny, sizeof(double *))))
    errmess("calloc(corrim)");
  for (j=0; j<ny; j++)
    if (!(corrim[j]=(double *)calloc((size_t)nx, sizeof(double))))
      errmess("calloc(corrim[j])");

  for (i=0; i<ny; i++)
    for (j=0; j<nx; j++)
      corrim[i][j]=0.0;

/* concolve and normalize to get covariance */
  for (j=0; j<ny; j++)
  {
    for (i=0; i<nx; i++)
    {
      if (im[j][i] > thresh)
      {
        norm=0;
        b = c = mean_im = 0.0;

        for (l=-par.mohw; l<=par.mohw; l++)
        {
          for(k=-par.mohw; k<=par.mohw; k++)
          {
            if ((j+l >= 0) && (j+l < ny) && (i+k >= 0) && (i+k < nx))
            {
              mean_im+=im[j+l][i+k];
              norm++;
            }
          }
        }

        mean_im/=norm;

        for (l=-par.mohw; l<=par.mohw; l++)
        {
          for (k=-par.mohw; k<=par.mohw; k++)
          {
            if ((j+l >= 0) && (j+l < ny) && (i+k >= 0) && (i+k < nx))
            {
              a0=im[j+l][i+k]-mean_im;
              c+=a0*a0;
              b+=(im[j+l][i+k]-mean_im)*(gauss[l][k]-mean_model);
            }
          }
        }

        if (c > 0.0)  corrim[j][i]=b/sqrt(a*c);
        else  if (par.verbose) printf("WARNING! covar: (%d, %d) c=0\n", i, j);
      }
    }
  }

  for (j=-par.mohw; j<=par.mohw; j++) free(gauss[j]-par.mohw);
  free(gauss-par.mohw);

  return(corrim);
}
