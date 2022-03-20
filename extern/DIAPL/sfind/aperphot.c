/*========================================================*/
/*                                                        */
/*  aperphot.c          version 1.6.2   2006.10.23        */
/*                                                        */
/*  Original source code:                                 */
/*  Copyright (C) 2006 by Przemek Wozniak                 */
/*  wozniak@lanl.gov                                      */
/*                                                        */
/*  Modifications:                                        */
/*  Copyright (C) 2006 by Wojtek Pych, CAMK PAN           */
/*  pych@camk.edu.pl                                      */
/*                                                        */
/*  Written for GNU project C and C++ Compiler            */
/*                                                        */
/* Given image im(nx,ny), centroid (cx,cy)                */
/* and background level calculates rough estimate         */
/* of the aperture photometry within radius rad.          */
/*                                                        */
/* sat_test= 0    OK                                      */
/* sat_test= 1    pixels above sat_level within aperture  */
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

float aperphot(int nx, int ny, float **im, int cx, int cy, float bglev,
                int *sat_test, PARAMS par)
{
        int   i, j, irad;
        float rr, rad, phot;

  *sat_test=0;
  irad=(int)par.aprad+1;
  rad=par.aprad*par.aprad;

  phot=0.0;

  for (j=-irad; j<=irad; j++)
  {
    if ((cy+j >= 0) && (cy+j < ny))
    {
      for (i=-irad; i<=irad; i++)
      {
        if ((cx+i >= 0) && (cx+i < nx))
        {
          rr=(float)i*i + (float)j*j;

          if (rr < rad)                       phot+=im[cy+j][cx+i]-bglev;
          if (im[cy+j][cx+i] > par.sat_level) *sat_test=1;
        }
      }
    }
  }

/* primitive magnitude zero point */
  if (phot > 0.0) return(38.0-2.5*log(phot));
  else            return(-999.999); // very low probability of a good result
}
