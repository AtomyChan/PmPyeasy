/*========================================================*/
/*                                                        */
/*  psf_core.c          version 1.6.2   2005.01.24        */
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

double psf_core(PSF_STRUCT *psf, double x, double y)
{
        int     m, n, ldeg, igauss, ngauss, icomp;
        double  f, x1, y1, a1, a2, rr, psf_pix;

  ngauss = psf->ngauss;
  ldeg   = psf->ndeg_local;

  x1 = psf->cos*x - psf->sin*y;
  y1 = psf->sin*x + psf->cos*y;

  rr = psf->ax*x1*x1 + psf->ay*y1*y1;

  psf_pix = 0.0;
  icomp = 0;

  for (igauss=0; igauss<ngauss; igauss++)
  {
    f = exp(rr);
    a1 = 1.0;
    for (m=0; m<=ldeg; m++)
    {
      a2 = 1.0;
      for (n=0; n<=ldeg-m; n++)
      {
        psf_pix += psf->local_vec[icomp++]*f*a1*a2;
        a2 *= y;
      }
      a1 *= x;
    }
    rr *= psf->sigma_inc*psf->sigma_inc;
  }

  return(psf_pix);
}
