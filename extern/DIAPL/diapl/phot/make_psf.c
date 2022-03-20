/*========================================================*/
/*                                                        */
/*  make_psf.c          version 1.6.2   2005.01.24        */
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

#include "defs.h"
#include "funcs.h"

void make_psf(PSF_STRUCT *psf, double xpsf, double ypsf, int cxpsf, int cypsf,
              double *psfim)
{
        int     cx, cy, psfhw, idx, size;
        double  x, y, dx, dy, rr, rad2;

  rad2  = psf->normrad;
  rad2 *= rad2;

  psfhw = psf->hw;
  size  = 2*psfhw + 1;

  dx = xpsf - (double)cxpsf;
  dy = ypsf - (double)cypsf;

  for (cx=-psfhw; cx<=psfhw; cx++)
  {
    for (cy=-psfhw; cy<=psfhw; cy++)
    {
      x  = (double)cx - dx;
      y  = (double)cy - dy;
      rr = x*x + y*y;

      idx = cx+psfhw + size*(cy+psfhw);

      if (rr <= rad2) psfim[idx] = psf_core(psf, x, y);
      else            psfim[idx] = 0.0;
    }
  }

  return;
}
