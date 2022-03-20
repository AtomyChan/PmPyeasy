/*========================================================*/
/*                                                        */
/*  base_func.c         version 1.3.0   2005.03.03        */
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

#include <stdio.h>
#include <math.h>

void base_func(double sigma, int m, int n, double *ker_x, double *ker_y,
               int kerhw, int do_norm)
{
        int     i, idx;
        double  sig, norm_x, norm_y ;

  sig = (double)(2*sigma*sigma);
  norm_x = norm_y = 0.0;

  for (i=-kerhw; i<=kerhw; i++)
  {
    idx = kerhw + i;

    ker_y[idx] = ker_x[idx]  = exp(-i*i/sig);

    if (m > 0)  ker_x[idx] *= pow(i, m);
    if (n > 0)  ker_y[idx] *= pow(i, n);

    norm_x += ker_x[idx];
    norm_y += ker_y[idx];
  }

  if (do_norm)
  {
    if ((norm_x != 0.0) && (norm_y != 0.0))
    {
      for (i=0; i<=2*kerhw; i++)
      {
        ker_x[i] /= norm_x;
        ker_y[i] /= norm_y;
      }
    }
    else
    {
      fprintf(stderr, "base_func: warning: normalization to 0.0 ignored!\n");
    }
  }

  return;
}
