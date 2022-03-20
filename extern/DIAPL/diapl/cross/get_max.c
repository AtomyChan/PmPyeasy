/*========================================================*/
/*                                                        */
/*  get_max.c           version 1.3.3   2005.01.28        */
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

void get_max(double *im1, double * im2, int nx, int ny,
             double *x_shift, double *y_shift)
{
        int   i, j, cx, cy;
        double marg[3], denom, dx, dy, max;

  for (j=0; j<ny/2; j++)
  {
    for (i=0; i<nx/2; i++)
    {
      im2[i+nx*j] = im1[2*(i+nx/2+nx*(j+ny/2))];
      im2[i+nx/2+nx*j] = im1[2*(i+nx*(j+ny/2))];
      im2[i+nx*(j+ny/2)] = im1[2*(i+nx/2+nx*j)];
      im2[i+nx/2+nx*(j+ny/2)] = im1[2*(i+nx*j)];
    }
  }

  max = -32000.0;
  cx = nx/2;
  cy = ny/2;

  for (j=1; j<ny-1; j++)
  {
    for (i=1; i<nx-1; i++)
    {
      if (im2[i+nx*j] > max)
      {
	max = im2[i+nx*j];
	cx = i;
	cy = j;
      }
    }
  }

/* marginal distributions */

  for (i=-1; i<=1; i++)
  {
    marg[i+1] = 0.0;

    for (j=-1; j<=1; j++)
    {
      marg[i+1] += im2[cx+i+nx*(cy+j)];
    }
  }

/* calculate position of the extremum of parabola through 3 points */

  denom = (marg[2] + marg[0] - 2.0*marg[1]);

  if(denom != 0.0)    dx = 0.5*(marg[0] - marg[2])/denom;
  else                dx = 0.0;

/* the same in y direction */

  for (j=-1; j<=1; j++)
  {
    marg[j+1] = 0.0;

    for (i=-1; i<=1; i++)
    {
      marg[j+1] += im2[cx+i+nx*(cy+j)];
    }
  }

  denom = (marg[2] + marg[0] - 2.0*marg[1]);

  if (denom != 0.0)    dy = 0.5*(marg[0] - marg[2])/denom;
  else                 dy = 0.0;

  if (fabs(dx) > 1.0) dx=0.0;
  if (fabs(dy) > 1.0) dy=0.0;

  *x_shift = dx + (double)cx;
  *y_shift = dy + (double)cy;

  return;
}
