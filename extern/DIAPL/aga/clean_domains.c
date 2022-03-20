/*========================================================*/
/*                                                        */
/*  clean_domains.c     version 1.3.0   2005.03.03        */
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

#include "defs.h"

int clean_domains(float *im, float *imref,
                  unsigned short *mask0, unsigned short *mask1,
                  double **vecs, DOM_TABS *domp, PARAMS par)
{
        int    i, j, m, n, idx, npix, xc, yc, idom;
        double **mat0, **mat1, *vec0;

  npix = 0;

  for (idom=0; idom<par.ndom; idom++)
  {
    mat0 = domp[idom].mat0;
    mat1 = domp[idom].mat1;
    vec0 = domp[idom].vec0;
    xc   = domp[idom].xc;
    yc   = domp[idom].yc;

    for (m=0; m<par.nvecs; m++)
    {
      vec0[m] = 0.0;
      for (n=0; n<=m; n++)  mat1[m][n] = mat0[m][n];
    }

    for (i=xc-par.domhw; i<=xc+par.domhw; i++)
    {
      for (j=yc-par.domhw; j<=yc+par.domhw; j++)
      {
        idx = i + par.nx*j;

        if (mask0[idx])
        {
          if (mask1[idx])
          {
            for (m=0; m<par.nvecs; m++)
              vec0[m] += im[idx]*vecs[m][idx]*imref[idx];
            npix++;
          }
          else
          {
            for (m=0; m<par.nvecs; m++)
              for (n=0; n<=m; n++)
                mat1[m][n] -= vecs[m][idx]*vecs[n][idx]*imref[idx];
          }
        }
      }
    }
  }

/* mask1 becomes logical product of mask0 and mask1 */

  for (idx=0; idx<par.nx*par.ny; idx++)
    if (!mask0[idx])  mask1[idx] = 0;

  return(npix);
}
