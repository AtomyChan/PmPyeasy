/*========================================================*/
/*                                                        */
/*  spatial_convolve.c  version 1.3.0   2005.03.03        */
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

void spatial_convolve(double *im, double *cnvim, double **vecs,
                      double **wxy, double *exp_vec, PARAMS par)
{
        int i, j, k, m, n, idx;

  for (i=par.kerhw; i<=par.nx-par.kerhw; i++)
  {
    for (j=par.kerhw; j<=par.ny-par.kerhw; j++)
    {
      idx = i + par.nx*j;

      cnvim[idx] = im[idx];

      for (k=0; k<par.nvecs; k++)
        cnvim[idx] -= (double)(exp_vec[k]*vecs[k][idx]);

      for (m=1; m<par.nwxy; m++)
      {
        for (n=par.nbkg+1; n<par.nvecs; n++)
        {
          cnvim[idx] -= (double)(exp_vec[k]*wxy[m][idx]*vecs[n][idx]);
          k++;
        }
      }
    }
  }

  return;
}
