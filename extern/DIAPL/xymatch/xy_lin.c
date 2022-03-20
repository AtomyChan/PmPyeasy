/*========================================================*/
/*                                                        */
/*  xy_lin.c            version 1.3.3   2005.02.10        */
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
/* Given (x1, y1) coordinates of nobj objects             */
/* and (x2, y2) coordinates of the same objects           */
/* on a different grid calculate first order coordinate   */
/* transformation between two frames. coeffx and coeffy   */
/* contain 3 coefficients each (in order x, y and shift). */
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

void xy_lin(float *x1, float *y1, float *x2, float *y2, int nobj,
            float *coeffx, float *coeffy)
{
        int   i;
        float sum_x2x2=0.0, sum_x2y2=0.0, sum_y2y2=0.0,
              sum_x2=0.0, sum_y2=0.0, sum_x1x2=0.0, sum_x1y2=0.0,
              sum_x1=0.0, sum_y1=0.0, sum_y1x2=0.0, sum_y1y2=0.0,
              a11, a12, a21, a22, b1, b2;

  for (i=0; i<nobj; i++)
  {
    sum_x2x2 += x2[i]*x2[i];
    sum_x2y2 += x2[i]*y2[i];
    sum_y2y2 += y2[i]*y2[i];

    sum_x1x2 += x1[i]*x2[i];
    sum_y1y2 += y1[i]*y2[i];
    sum_x1y2 += x1[i]*y2[i];
    sum_y1x2 += y1[i]*x2[i];

    sum_x1 += x1[i];
    sum_x2 += x2[i];
    sum_y1 += y1[i];
    sum_y2 += y2[i];
  }

  a11 = sum_x2x2 - sum_x2*sum_x2/nobj;
  a12 = sum_x2y2 - sum_x2*sum_y2/nobj;
  a21 = sum_x2y2 - sum_x2*sum_y2/nobj;
  a22 = sum_y2y2 - sum_y2*sum_y2/nobj;

  b1  = sum_x1x2 - sum_x1*sum_x2/nobj;
  b2  = sum_x1y2 - sum_x1*sum_y2/nobj;

  coeffx[0] = ( b1/a12 - b2/a22 ) / ( a11/a12 - a21/a22 );
  coeffx[1] = ( b1/a12 - coeffx[0]*a11/a12 );
  coeffx[2] = ( sum_x1 - coeffx[0]*sum_x2 - coeffx[1]*sum_y2 )/nobj;

  b1  = sum_y1x2 - sum_y1*sum_x2/nobj;
  b2  = sum_y1y2 - sum_y1*sum_y2/nobj;

  coeffy[0] = ( b1/a12 - b2/a22 ) / ( a11/a12 - a21/a22 );
  coeffy[1] = ( b1/a12 - coeffy[0]*a11/a12 );
  coeffy[2] = ( sum_y1 - coeffy[0]*sum_x2 - coeffy[1]*sum_y2 )/nobj;

  return;
}
