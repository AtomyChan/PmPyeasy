/*========================================================*/
/*                                                        */
/*  poly.c              version 1.2.5   2005.04.07        */
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
/* Evaluates polynomial of order ndeg in x and y          */
/* with coefficients stored in coeff.                     */
/*                                                        */
/*========================================================*/

#include <stdlib.h>

#include "errmess.h"

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

float poly(float x, float y, int ndeg, double *coeff)
{
        int   j, k, njk, nterm;
        float pol, *f;

  pol   = coeff[0];
  nterm = (ndeg+1)*(ndeg+2)/2;

  if (!(f=(float *)calloc(nterm, sizeof(float)))) errmess("calloc(f)");

  f[0] = 1.0;
  njk  = 1;

  for (j=1; j<=ndeg; j++)
  {
    for (k=0; k<j; k++)
    {
      f[njk] = f[njk-j]*y;
      pol   += coeff[njk]*f[njk];
      ++njk;
    }

    f[njk] = f[njk-j-1]*x;
    pol   += coeff[njk]*f[njk];
    ++njk;
  }

  free(f);

  return(pol);
}
