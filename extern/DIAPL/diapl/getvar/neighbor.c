/*========================================================*/
/*                                                        */
/*  neighbor.c          version 1.2.3   2005.01.24        */
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

int neighbor(double *im, int cx, int cy, PAR_STRUCT *par, double bg)
{
        int   xlo, xhi, ylo, yhi, i, j, k, l, is_max, nx, ny;
        double pix, pix0, rad, nsig, gain;

  nx = par->nx0;
  ny = par->ny0;

  gain = par->gain;

  xlo = cx-par->isohw;
  xhi = cx+par->isohw;
  ylo = cy-par->isohw;
  yhi = cy+par->isohw;

  if (xlo < par->smhw)  xlo=par->smhw;
  if (ylo < par->smhw)  ylo=par->smhw;
  if (xhi >= nx-par->smhw)  xhi=nx-par->smhw-1;
  if (yhi >= ny-par->smhw)  yhi=nx-par->smhw-1;

  pix0 = im[cx+nx*cy];

  for (i=xlo; i<=xhi; i++)
  {
    for (j=ylo; j<=yhi; j++)
    {
      pix = im[i+nx*j];

      is_max = 1;

      for (k=i-par->smhw; k<=i+par->smhw; k++)
        for (l=j-par->smhw; l<=j+par->smhw; l++)
          if (pix < im[k+nx*l])
            is_max = 0;

      nsig = (pix-bg)/sqrt(pix/gain);
      rad  = sqrt((i-cx)*(i-cx) + (j-cy)*(j-cy));

      if (is_max && (rad >= 1.0) && (nsig >= par->nsig_var1))
        if (pix > (par->iso_off + par->iso_slo*rad)*pix0)
          return(1);
    }
  }

  return(0);
}
