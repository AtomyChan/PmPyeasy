/*========================================================*/
/*                                                        */
/*  get_fwhm.c          version 1.2.3   2005.01.24        */
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

#include "defs.h"

double get_fwhm(double *psfim, double x, double y, int cx, int cy,
                PAR_STRUCT *par, double *ratio)
{
        int     i, j, kx, ky, irad, psfn, psfhw;
        double  rr, rad2, psf_pix, sxx, sxy, syy, sum, delta,
                sigma_x, sigma_y, dx, dy, fwhm;

  rad2 = par->fitrad*par->fitrad;

  irad  = (int)par->fitrad + 2;
  psfhw = par->psfhw;
  psfn  = par->psfn;

  sxx = sxy = syy = sum = 0.0;

  for (i=cx-irad; i<=cx+irad; i++)
  {
    for( j=cy-irad; j<=cy+irad; j++)
    {
      kx = i - cx;
      ky = j - cy;

      dx = i - x;
      dy = j - y;

      rr = dx*dx + dy*dy;

      psf_pix = psfim[psfhw+kx+psfn*(psfhw+ky)];

      if (rr <= rad2)
      {
      	sxx += psf_pix*dx*dx;
        sxy += psf_pix*dx*dy;
        syy += psf_pix*dy*dy;
        sum += psf_pix;
      }
    }
  }

  if (sum != 0.0)
  {
    sxx /= sum;
    sxy /= sum;
    syy /= sum;
  }

  delta = sxx*sxx - 2.0*sxx*syy + syy*syy + 4.0*sxy*sxy;

  sigma_x = (sxx + syy + sqrt(delta))*0.5;
  sigma_y = (sxx + syy - sqrt(delta))*0.5;

  if ((sigma_y > 0.0) && (sigma_x > 0.0))
  {
    if (sigma_y > sigma_x)  *ratio = sqrt(sigma_x/sigma_y);
    else                    *ratio = sqrt(sigma_y/sigma_x);

    fwhm = sqrt(4.0*log(2.0)*(sigma_x + sigma_y));

//    if (fwhm > 5.0) fwhm = 5.0;
  }
  else
  {
    *ratio = 0.0;
    fwhm = 99.9;
  }

  return (fwhm);
}
