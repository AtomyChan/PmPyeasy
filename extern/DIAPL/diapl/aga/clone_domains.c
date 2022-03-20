/*========================================================*/
/*                                                        */
/*  clone_domains.c     version 1.3.0   2005.03.03        */
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

void clone_domains(DOM_TABS *domp, int ndom, int nvecs)
{
        int    m, n, idom;
        double **mat, **mat1, *vec0, *vec;

  for (idom=0; idom<ndom; idom++)
  {
    mat  = domp[idom].mat;
    vec  = domp[idom].vec;
    mat1 = domp[idom].mat1;
    vec0 = domp[idom].vec0;

    for (m=0; m<nvecs; m++)
    {
      vec[m] = vec0[m];
      for (n=0; n<=m; n++)  mat[m][n] = mat1[m][n];
    }

    for (m=0; m<nvecs; m++)
      for(n=0; n<m; n++) mat[n][m] = mat[m][n];
  }

  return;
}
