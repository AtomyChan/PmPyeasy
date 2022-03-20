/*========================================================*/
/*                                                        */
/*  init_psf.c          version 1.6.2   2005.01.24        */
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

void init_psf(PSF_STRUCT *psf, double xpsf, double ypsf)
{
        int     m, n, icomp, ncomp, itot, sdeg, ldeg;
        double  a1, a2, *vec, *local_vec, x_orig, y_orig;

  x_orig = psf->x_orig;
  y_orig = psf->y_orig;

  sdeg   = psf->ndeg_spat;
  ldeg   = psf->ndeg_local;
  ncomp  = psf->ngauss*(ldeg+1)*(ldeg+2)/2;

  local_vec = psf->local_vec;
  vec = psf->vec;

  for (icomp=0; icomp<ncomp; icomp++)
    local_vec[icomp] = 0.0;

  itot = 0;
  a1 = 1.0;

  for (m=0; m<=sdeg; m++)
  {
    a2 = 1.0;
    for (n=0; n<=sdeg-m; n++)
    {
      for (icomp=0; icomp<ncomp; icomp++)
        local_vec[icomp] += vec[itot++]*a1*a2;

      a2 *= ypsf - y_orig;
    }
    a1 *= xpsf - x_orig;
  }

  return;
}
