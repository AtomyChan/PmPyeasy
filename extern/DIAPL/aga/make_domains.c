/*========================================================*/
/*                                                        */
/*  make_domains.c      version 1.3.0   2005.03.03        */
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

void make_domains(float *imref, unsigned short *mask0, double **vecs,
                  DOM_TABS *domp, PARAMS par)
{
        int i, j, m, n, idx, xc, yc, idom;
        double **mat0;

/* avoid divisions in future calculations */
  for (idx=0; idx<par.nx*par.ny; idx++)
    if (mask0[idx]) imref[idx] = 1.0/imref[idx];

  for (idom=0; idom<par.ndom; idom++)
  {
    mat0 = domp[idom].mat0;
    xc   = domp[idom].xc;
    yc   = domp[idom].yc;

    for (m=0; m<par.nvecs; m++)
      for (n=0; n<=m; n++)
        mat0[m][n] = 0.0;

    for (i=xc-par.domhw; i<=xc+par.domhw; i++)
    {
      for (j=yc-par.domhw; j<=yc+par.domhw; j++)
      {
        idx = i + par.nx*j;

        if (mask0[idx])
        {
          for (m=0; m<par.nvecs; m++)
            for (n=0; n<=m; n++)
              mat0[m][n] += vecs[m][idx]*vecs[n][idx]*imref[idx];
        }
      }
    }
  }

  return;
}
