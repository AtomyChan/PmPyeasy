/*========================================================*/
/*                                                        */
/*  is_isolated.c       version 1.3.1   2005.01.28        */
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

int is_isolated(IMAGE *im, OBJECT *obj, FITPAR par)
{
        int   cx, cy, xlo, xhi, ylo, yhi, i, j, k, l, is_max, nx, ny;
        double pix, pix0, rad, bg, nsig;

  nx = im->nx;
  ny = im->ny;

  cx  = obj->x;
  cy  = obj->y;
  bg  = obj->bg;

  xlo = cx-par.isohw;
  xhi = cx+par.isohw;
  ylo = cy-par.isohw;
  yhi = cy+par.isohw;

  if (xlo<par.maxhw)  xlo=par.maxhw;
  if (ylo<par.maxhw)  ylo=par.maxhw;
  if (xhi>=nx-par.maxhw)  xhi=nx-par.maxhw-1;
  if (yhi>=ny-par.maxhw)  yhi=nx-par.maxhw-1;

  pix0 = im->pix[cx+nx*cy];

  for (i=xlo; i<=xhi; i++)
  {
    for (j=ylo; j<=yhi; j++)
    {
      pix = im->pix[i+nx*j];

      is_max = 1;

      for (k=i-par.maxhw; k<=i+par.maxhw; k++)
        for (l=j-par.maxhw; l<=j+par.maxhw; l++)
          if (pix < im->pix[k+nx*l])  is_max = 0;

      nsig = (pix-bg)/sqrt(pix/im->gain);
      rad  = sqrt((i-cx)*(i-cx) + (j-cy)*(j-cy));

      if (is_max && (rad >= 1.0) && (nsig >= par.nsig_detect))
        if (pix > (par.iso_off + par.iso_slo*rad)*pix0)
          return(0);
    }
  }

  return(1);
}
