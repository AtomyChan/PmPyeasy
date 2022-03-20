/*========================================================*/
/*                                                        */
/*  fitn.c              version 1.1.4   2005.03.09        */
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
/* Fits z as a polynomial of order ndeg in x and y.       */
/* All coordinate lists are ndata long.                   */
/* v contains coefficients.                               */
/*                                                        */
/* uses numerical recipes for matrix inversion            */
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
#include <stdlib.h>

#include "errmess.h"

#include "ludcmp.h"
#include "lubksb.h"

void fit(double *x, double *y, double *z, int ndeg, int ndata, double *v)
{
        int     i, j, k, njk, nterm, *index;
        double   *f;
        double  **a, d;

  nterm = (ndeg+1)*(ndeg+2)/2;

  if (!(index=(int *)malloc(nterm*sizeof(int)))) errmess("malloc(index)");

/* beware of stupid 1,n indexing of numerical recipes ! */

  if (!(a=(double **)malloc((nterm+1)*sizeof(double *)))) errmess("malloc(a)");
  for (i=0; i<=nterm; i++)
    if (!(a[i]=(double *)malloc((nterm+1)*sizeof(double))))
      errmess("malloc(a[i])");

  if (!(f=(double *)malloc(nterm*sizeof(double)))) errmess("malloc(f)");

  for (j=0; j<nterm; j++)
  {
    v[j] = 0.0;

    for (k=0; k<nterm; k++) a[j+1][k+1] = 0.0;
  }

/* build moment table */

  for (i=0; i<ndata; i++)
  {
    f[0] = 1.0;
    njk  = 1;

    for (j=1; j<=ndeg; j++)
    {
      for (k=0; k<j; k++)
      {
        f[njk] = f[njk-j]*y[i];
        ++njk;
      }

      f[njk] = f[njk-j-1]*x[i];
      ++njk;
    }

    for (j=0; j<nterm; j++)
    {
      v[j] += z[i]*f[j];

      for (k=0; k<=j; k++)  a[j+1][k+1] += f[j]*f[k];
    }
  }

  for (j=1; j<=nterm; j++)
    for(k=1; k<j; k++)
      a[k][j] = a[j][k];

/* solve linear system of least square equations by matrix inversion */
/* beware of stupid 1,n indexing of numerical recipes !              */

  ludcmp(a, nterm, index-1, &d);
  lubksb(a, nterm, index-1, v-1);

  free(f);
  free(index);

  for (i=0; i<=nterm; i++) free(a[i]);
  free(a);

  return;
}
