/*========================================================*/
/*                                                        */
/*  center.c            version 1.3.1   2005.01.28        */
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
/* Given im(nx, ny) and pixel (cx, cy) which is a local   */
/* maximum calculates centroid (x,y).                     */
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
#include <stdio.h>
#include <stdlib.h>

#include "defs.h"

void center(IMAGE *im, OBJECT *obj)
{
        int   cx, cy, i, j, nx, ny;
        double marg[3], dx, dy, denom;

  nx = im->nx;
  ny = im->ny;

  cx = obj->x;
  cy = obj->y;

  if ((cx < 1) || (cx >= nx-1) || (cy < 1) || (cy >= ny-1))
  {
    fprintf(stderr, "centroid: pixel out of bounds!\n");
    exit(1);
  }

/* take projection of 3x3 pix neighborhood */

  for (i=-1; i<=1; i++)
  {
    marg[i+1]=0.0;

    for (j=-1; j<=1; j++) marg[i+1] += im->pix[cx+i+nx*(cy+j)];
  }

/* calculate position of the extremum of parabola through 3 points */

  denom = (marg[2] + marg[0] - 2.0*marg[1]);

  if (denom != 0.0) dx=0.5*(marg[0] - marg[2])/denom;
  else              dx=0.0;

/* the same in y direction */

  for (j=-1; j<=1; j++)
  {
    marg[j+1]=0.0;

    for (i=-1; i<=1; i++) marg[j+1] += im->pix[cx+i+nx*(cy+j)];
  }

  denom = (marg[2] + marg[0] - 2.0*marg[1]);

  if (denom != 0.0) dy=0.5*(marg[0] - marg[2])/denom;
  else              dy=0.0;

  if (fabs(dx) > 1.0) dx=0.0;
  if (fabs(dy) > 1.0) dy=0.0;

  obj->xfit = dx + (double)cx;
  obj->yfit = dy + (double)cy;

  return;
}
