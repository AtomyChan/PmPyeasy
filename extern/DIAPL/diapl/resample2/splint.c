/*========================================================*/
/*                                                        */
/*  splint.c            version 1.2.4   2005.02.16        */
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
/* Spline evaluation for entire vector in one go,         */
/* rewritten from numerical recipes (uses 0,n-1 indexing) */
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

void splint(double *xa, double *ya, double *y2a, int n, double *x, double *y,
            RINGPAR par)
{
        int   i, k;
        double h, b[512], a[512], Dy;

  k = 0;
  for (i=0; i<n; i++)
  {
    while (x[i] > xa[k]) k++;

    if (k < 1)   y[i] = 1111.1111;
    else if   (k > n-1)                          y[i] = 1234.0;
    else if ((ya[k] < par.min) || (ya[k-1] < par.min) ||
             (ya[k] > par.sat) || (ya[k-1] > par.sat))  y[i] = 2345.0;
    else
    {
      h = xa[k] - xa[k-1];
      if (h < 0.0)
        fprintf(stderr, "splint: warning: vector XA should be sorted\n");
      a[i] = (xa[k] - x[i])/h;
      b[i] = (x[i] - xa[k-1])/h;
      y[i] = a[i]*ya[k-1] + b[i]*ya[k];
      Dy   = ( (a[i]*a[i]*a[i]-a[i])*y2a[k-1] + (b[i]*b[i]*b[i]-b[i])*y2a[k] )*(h*h)/6.0;

      if (par.fix)
      {
        if (Dy*Dy < par.nsig*par.nsig*ya[i]/par.gain) y[i] += Dy;
      }
      else                                            y[i] += Dy;
    }
  }
//    printf("to check the splint function\n");
//    for (i=0;i<5;i++) printf("a[%d]=%g   b[%d]=%g\n",i,a[i],i,b[i]);
  return;
}
