/*========================================================*/
/*                                                        */
/*  get_gaussian.c      version 1.3.1   2005.01.28        */
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

#include <stdio.h>
#include <math.h>

#include "defs.h"
#include "funcs.h"

void get_gaussian(PSF_STRUCT *psf)
{
        int     cx, cy, m, n, psfhw, icomp, igauss, ngauss, ncomp, ldeg;
        double  a1, a2, *vec, psf_pix, x1, y1, x2, y2, rr, sum, sxx, syy, sxy,
                theta_psf=0.0, delta, sigma2_mscale, sigma2_inc, f;

  psfhw   = psf->hw;
  vec     = psf->vec;
  ldeg    = psf->ndeg_local;
  ngauss  = psf->ngauss;

  sigma2_inc    = psf->sigma_inc*psf->sigma_inc;
  sigma2_mscale = psf->sigma_mscale*psf->sigma_mscale;

  ncomp  = ngauss*(ldeg+1)*(ldeg+2)/2;

  sum = sxx = syy = sxy = 0.0;

  for (cx=-psfhw; cx<=psfhw; cx++)
  {
    for (cy=-psfhw; cy<=psfhw; cy++)
    {
      x1 = (double)cx;
      y1 = (double)cy;

      x2 = psf->cos*x1 - psf->sin*y1;
      y2 = psf->sin*x1 + psf->cos*y1;

      rr = psf->ax*x2*x2 + psf->ay*y2*y2;

      icomp = 0;
      psf_pix = 0.0;

      for (igauss=0; igauss<ngauss; igauss++)
      {
        f = exp(rr);
        a1 = 1.0;
        for (m=0; m<=ldeg; m++)
        {
          a2 = 1.0;
          for (n=0; n<=ldeg-m; n++)
          {
            psf_pix += vec[icomp++]*f*a1*a2;
            a2 *= y1;
          }
          a1 *= x1;
        }
        rr *= sigma2_inc;
      }

      sxx += psf_pix*x1*x1;
      sxy += psf_pix*x1*y1;
      syy += psf_pix*y1*y1;
      sum += psf_pix;
    }
  }

  if (sum != 0.0)
  {
    sxx /= sum;
    sxy /= sum;
    syy /= sum;
  }

  delta = sxx*sxx - 2.0*sxx*syy + syy*syy + 4.0*sxy*sxy;

  psf->ax = (sxx + syy + sqrt(delta))*0.5;
  psf->ay = (sxx + syy - sqrt(delta))*0.5;

  if (sxy != 0.0) theta_psf = atan((2.0*psf->ax - sxx)/sxy);

  psf->cos = cos(theta_psf);
  psf->sin = sin(theta_psf);
  psf->ax  = -1.0/(2.0*psf->ax)*sigma2_mscale;
  psf->ay  = -1.0/(2.0*psf->ay)*sigma2_mscale;

  return;
}
