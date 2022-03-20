/*========================================================*/
/*                                                        */
/*  aperphot.c          version 1.2.3   2005.01.24        */
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
/* Given image im(nx, ny), centroid (cx, cy)              */
/* and background level, calculates rough estimate        */
/* of the aperture photometry within radius rad.          */
/* sat_test = 0 in case there are pixels above sat_level  */
/*              within aperture.                          */
/* sat_test = 1 otherwise                                 */
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

double aperphot(double *im, int nx, int ny, int cx, int cy, double rad,
               double bglev, int *sat_test, double sat_level)
{
        int   i, j, irad, idx;
        double rr, phot;

  irad = rad + 1;

  phot = 0.0;

  for (j=-irad; j<=irad; j++)
  {
    for (i=-irad; i<=irad; i++)
    {
      rr  = i*i + j*j;
      idx = cx+i+nx*(cy+j);

      if (rr < rad*rad) phot  += im[idx] - bglev;
      if (im[idx] > sat_level) *sat_test=0;
    }
  }

  return(phot);
}
