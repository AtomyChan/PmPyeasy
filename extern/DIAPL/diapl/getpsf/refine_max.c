/*========================================================*/
/*                                                        */
/*  refine_max.c        version 1.3.1   2005.01.28        */
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

#include "defs.h"

void refine_max(IMAGE *im, OBJECT *obj)
{
        int    k, xc, yc, nx;
        double *pix, a0, a1, a2, ym[3], max;

  pix = im->pix;

  xc = obj->x;
  yc = obj->y;

  nx = im->nx;

  for (k=-1; k<=1; k++)
  {
    a0 =  pix[xc  +nx*(yc+k)];
    a1 = (pix[xc+1+nx*(yc+k)] - pix[xc-1+nx*(yc+k)])/2.0;
    a2 = (pix[xc+1+nx*(yc+k)] + pix[xc-1+nx*(yc+k)])/2.0 - a0;

    if (a2!=0.0)  ym[k+1] = a0 - a1*a1/(4.0*a2);
    else          ym[k+1] = a0;
  }

  a0 =  ym[1];
  a1 = (ym[2] - ym[0])/2.0;
  a2 = (ym[2] + ym[0])/2.0 - a0;

  if (a2!=0.0)  max = a0 - a1*a1/(4.0*a2);
  else          max = a0;

  obj->max = max;

  return;
}
