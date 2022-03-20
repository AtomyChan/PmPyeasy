/*========================================================*/
/*                                                        */
/*  get_phot.c          version 1.6.2   2005.01.24        */
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

void get_phot(float *difim, float *im, float *refim, double x, double y,
              int cx, int cy, double *psfim, STAR *obj, PAR_STRUCT *par)
{
        int     i, j, kx, ky, irad, count, nbad, psfn, psfhw, nx0, idx;
        double  rr, norm, p_sig2, a_sig2, p_phot, a_phot,
                a, b, c, sum_a, sum_b, sum_c, p_sum2, rad2,
                bg, sat, min, gain, pix, pix_psf, weight, pix_ref;

  rad2 = par->fitrad*par->fitrad;
  gain = par->gain;
  sat  = par->sat_level;
  min  = par->min_level;
  bg   = obj->bg;

  irad  = (int)par->fitrad + 2;
  psfhw = par->psfhw;
  psfn  = par->psfn;
  nx0   = par->nx0;

  count = nbad   = 0;
  sum_a = sum_b = sum_c = 0.0;
  p_sum2 = p_phot = p_sig2 = 0.0;
  a_phot = a_sig2 = norm = 0.0;

  for (i=0; i<psfn*psfn; i++) norm += psfim[i];

  for (i=cx-irad; i<=cx+irad; i++)
  {
    for (j=cy-irad; j<=cy+irad; j++)
    {
      kx = i - cx;
      ky = j - cy;

      rr = (i-x)*(i-x) + (j-y)*(j-y);

      idx = i + nx0*j;

      pix_psf = psfim[psfhw+kx+psfn*(psfhw+ky)];
      pix     = difim[idx];
      weight  =    im[idx];
      pix_ref = refim[idx];

      if (rr <= rad2)
      {
        if ((pix != par->err) && (weight < sat) && (pix_ref < sat) &&
            (weight >= min) && (pix_ref >= min))
        {

          a = pix_psf*(pix - bg);
          b = pix_psf*pix_psf;
          c = (pix - bg)*(pix - bg);

          sum_a += a;
          sum_b += b;
          sum_c += c;

          a_phot += pix - bg;
          a_sig2 += weight;
          p_phot += a/(weight*norm);
          p_sig2 += b/(weight*norm*norm);
          p_sum2 += c/weight;

          count++;
        }
        else  nbad++;
      }
    }
  }

  if ((p_sig2 != 0.0) && (count > 0))
  {
    obj->p_flux = (float)(p_phot/p_sig2);
    obj->p_err  = (float)(1.0/sqrt(p_sig2*gain));
    obj->chi2_n = (float)((p_sum2 - p_phot*p_phot/p_sig2)*gain/count);
    obj->corr   = (float)(sum_a/sqrt(sum_b*sum_c));

    obj->a_flux = (float)(a_phot);
    obj->a_err  = (float)(sqrt(a_sig2/gain));
    obj->nbad   = nbad;
  }
  else
  {
    obj->p_flux = par->bad_value;
    obj->p_err  = par->bad_value;
    obj->chi2_n = par->bad_value;
    obj->corr   = par->bad_value;

    obj->a_flux = par->bad_value;
    obj->a_err  = par->bad_value;
    obj->nbad   = -1;
  }

  return;
}
/*** END ***/
