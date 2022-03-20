/*========================================================*/
/*                                                        */
/*  is_peak.c           version 1.3.2   2005.04.18        */
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

int is_peak(IMAGE *im, int xc, int yc, int hw, double contrast)
{
        int   i, j;
        double max, min, pix;

  max = im->pix[xc+im->nx*yc];
  min = max;

  for (i=xc-hw; i<=xc+hw; i++)
  {
    for (j=yc-hw; j<=yc+hw; j++)
    {
      pix = im->pix[i+im->nx*j];

      if ((pix > max) || (pix >= im->sat_level) || (pix < im->min_level))
	return(0);

      if (pix < min) min = pix;
    }
  }

  if ((max-min)*(max-min)*im->gain < contrast*contrast*max) return(0);
  else                                                      return(1);
}
