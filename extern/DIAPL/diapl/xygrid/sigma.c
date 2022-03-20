/*========================================================*/
/*                                                        */
/*  sigma.c             version 1.1.4   2005.03.09        */
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
/* Calculates std deviation around the ndeg order         */
/* polynomial fit to z in x and y.                        */
/* Coefficients stored in coeff and all lists ndata long. */
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

#include "poly.h"

double sigma(double *x, double *y, double *z, int ndeg, int ndata, double *coeff)
{
        int   i;
        double sig, sg;

  sig = 0.0;

  for (i=0; i<ndata; i++)
  {
    sg   = z[i] - poly(x[i], y[i], ndeg, coeff);
    sg  *= sg;
    sig += sg;
  }

  sig = sqrt(sig/(double)(ndata-1));

  return sig;
}
