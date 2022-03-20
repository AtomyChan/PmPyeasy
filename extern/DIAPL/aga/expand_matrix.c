/*========================================================*/
/*                                                        */
/*  expand_matrix.c     version 1.3.0   2005.03.03        */
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

void expand_matrix(double **exp_mat, double *exp_vec, double **wxy,
                   DOM_TABS *domp, PARAMS par)
{
        int     m, n, ki, kj, li, lj, xc, yc, idom, idx;
        double  **mat, *vec;

  for (m=0; m<par.ntot; m++)
  {
    exp_vec[m] = 0.0;
    for (n=0; n<=m; n++)  exp_mat[m+1][n+1] = 0.0;
  }

  for (idom=0; idom<par.ndom; idom++)
  {
    if (domp[idom].reject != RJCT_NONE) continue;

    vec = domp[idom].vec;
    mat = domp[idom].mat;
    xc  = domp[idom].xc;
    yc  = domp[idom].yc;

    idx = xc + par.nx*yc;

    for (m=0; m<par.nvecs; m++)
    {
      exp_vec[m] += vec[m];
      for (n=0; n<=m; n++)  exp_mat[m+1][n+1] += mat[m][n];
    }

    m = par.nvecs;
    for (ki=1; ki<par.nwxy; ki++)
    {
      for (li=par.nbkg+1; li<par.nvecs; li++)
      {
        exp_vec[m] += wxy[ki][idx]*vec[li];

        for (n=0; n<par.nvecs; n++)
          exp_mat[m+1][n+1] += wxy[ki][idx]*mat[li][n];

        m++;
      }
    }

    m = par.nvecs;
    for (ki=1; ki<par.nwxy; ki++)
    {
      for (li=par.nbkg+1; li<par.nvecs; li++)
      {
        n = par.nvecs;
        for (kj=1; kj<par.nwxy && n<=m; kj++)
        {
          for (lj=par.nbkg+1; lj<par.nvecs && n<=m; lj++)
          {
            exp_mat[m+1][n+1] += wxy[ki][idx]*wxy[kj][idx]*mat[li][lj];
            n++;
          }
        }
        m++;
      }
    }
  }

  for (m=0; m<par.ntot; m++)
    for (n=0; n<m; n++)
      exp_mat[n+1][m+1] = exp_mat[m+1][n+1];

  return;
}
