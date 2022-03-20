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

void splint(float *xa, float *ya, float *y2a, int n, float *x, float *y,
            RINGPAR par)
{
        int   i, k;
        float h, b, a, Dy;

  k = 0;
  for (i=0; i<n; i++)
  {
    while (x[i] > xa[k]) k++;

    if ((k < 1) || (k > n-1))                           y[i] = 0.0;
    else if ((ya[k] < par.min) || (ya[k-1] < par.min) ||
             (ya[k] > par.sat) || (ya[k-1] > par.sat))  y[i] = 0.0;
    else
    {
      h = xa[k] - xa[k-1];
      if (h < 0.0)
        fprintf(stderr, "splint: warning: vector XA should be sorted\n");
      a = (xa[k] - x[i])/h;
      b = (x[i] - xa[k-1])/h;
      y[i] = a*ya[k-1] + b*ya[k];
      Dy   = ( (a*a*a-a)*y2a[k-1] + (b*b*b-b)*y2a[k] )*(h*h)/6.0;

      if (par.fix)
      {
        if (Dy*Dy < par.nsig*par.nsig*ya[i]/par.gain) y[i] += Dy;
      }
      else                                            y[i] += Dy;
    }
  }

  return;
}
