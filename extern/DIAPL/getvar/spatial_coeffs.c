/*========================================================*/
/*                                                        */
/*  spatial_coeffs.c    version 1.2.3   2005.01.24        */
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

void spatial_coeffs(KER_STRUCT *ker, double x, double y, double *wxy)
{
        int     m, n, hnx, hny, nwxy, wdeg;
        double  px, py, ax, ay, dx, dy;

  wdeg = ker->wdeg;
  hnx  = ker->nx/2;
  hny  = ker->ny/2;

  dx = 1.0/(hnx - ker->hw);
  dy = 1.0/(hny - ker->hw);

  px = (x - hnx)*dx;
  py = (y - hny)*dy;

  nwxy = 0;
  for (m=0, ax=1.0; m<=wdeg; m++, ax*=px)
    for (n=0, ay=1.0; n<=wdeg-m; n++, ay*=py)
      wxy[nwxy++] = ax*ay;

  if (nwxy != ker->nwxy)
  {
    fprintf(stderr, "spacial_coeffs: warning: wrong number of ");
    fprintf(stderr, "spatial vectors - shouldn't happen!\n");
    fflush(stderr);
  }

  return;
}
