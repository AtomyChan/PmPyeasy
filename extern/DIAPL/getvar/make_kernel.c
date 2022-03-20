/*========================================================*/
/*                                                        */
/*  make_kernel.c       version 1.2.3   2005.01.24        */
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

void make_kernel(KER_STRUCT *ker, double *wxy, double *kerim, int sector)
{
        int     k, m, n, idx, nvecs, nwxy, kern;
        double  **vecs, *ker_vec;

  ker_vec = ker->vec[sector];
  vecs    = ker->vecs;

  nvecs = ker->nvecs;
  nwxy  = ker->nwxy;
  kern  = 2*ker->hw + 1;

  for (idx=0; idx<kern*kern; idx++)
  {
    kerim[idx] = 0.0;

    for (k=0; k<nvecs; k++) kerim[idx] += ker_vec[k]*vecs[k][idx];

    for (m=1; m<nwxy; m++)
      for (n=1; n<nvecs; n++)
        kerim[idx] += ker_vec[k++]*wxy[m]*vecs[n][idx];
  }

  return;
}
