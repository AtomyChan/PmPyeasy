/*========================================================*/
/*                                                        */
/*  spline.c            version 1.2.4   2005.02.16        */
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
/* Spline contruction rewritten from numerical recipes    */
/* to use 0,n-1 indexing                                  */
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

#include <stdlib.h>

#include "errmess.h"

void spline(float *x, float *y, int n, float yp1, float ypn, float *y2)
{
        int   i, k;
        float p, qn, sig, un, *u;

  if (!(u=(float *)malloc((n-1)*sizeof(float)))) errmess("spline: malloc(u)");

  if (yp1 > 0.99e30)
    y2[0] = u[0] = 0.0;
  else
  {
    y2[0] = -0.5;
    u[0] = (3.0/(x[1]-x[0]))*((y[1]-y[0])/(x[1]-x[0])-yp1);
  }

  for (i=1; i<n-1; i++)
  {
    sig = (x[i]-x[i-1])/(x[i+1]-x[i-1]);
    p = sig*y2[i-1]+2.0;

    y2[i] = (sig-1.0)/p;

    u[i] = (y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
    u[i] = (6.0*u[i]/(x[i+1]-x[i-1]) - sig*u[i-1])/p;
  }

  if (ypn > 0.99e30)
    qn = un = 0.0;
  else
  {
    qn = 0.5;
    un = (3.0/(x[n-1]-x[n-2]))*(ypn-(y[n-1]-y[n-2])/(x[n-1]-x[n-2]));
  }

  y2[n-1] = (un-qn*u[n-2])/(qn*y2[n-2]+1.0);

  for (k=n-2; k>=0; k--)  y2[k] = y2[k]*y2[k+1]+u[k];

  free(u);

  return;
}
