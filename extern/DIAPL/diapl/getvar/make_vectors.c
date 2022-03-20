/*========================================================*/
/*                                                        */
/*  make_vectors.c      version 1.2.3   2005.01.24        */
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
#include <stdlib.h>

#include "defs.h"
#include "funcs.h"
#include "errmess.h"

void make_vectors(KER_STRUCT *ker)
{
        int     *deg, nvecs, mode,
                i, j, m, n, k, kern, kerhw, idx, ncomp;
        double   *sig;
        double  *ker_x, *ker_y, **vecs;

  sig  = ker->sig;
  deg  = ker->deg;

  vecs  = ker->vecs;
  ncomp = ker->ncomp;
  kerhw = ker->hw;
  kern  = 2*kerhw + 1;

  if (!(ker_x=(double *)malloc(kern*sizeof(double)))) errmess("malloc(ker_x)");
  if (!(ker_y=(double *)malloc(kern*sizeof(double)))) errmess("malloc(ker_y)");

  nvecs = 0;

  for (k=0; k<ncomp; k++)
  {
    for (m=0; m<=deg[k]; m++)
    {
      for (n=0; n<=deg[k]-m; n++)
      {
        mode = ((m/2)*2-m+1)*((n/2)*2-n+1);

        base_func(sig[k], m, n, ker_x, ker_y, kerhw, mode);

/* the first vector is normalized to 1.0   */
/* the rest must have zero sums to enforce */
/* flux conservation                       */

        for (i=0; i<kern; i++)
          for (j=0; j<kern; j++)
            vecs[nvecs][i+kern*j] = ker_x[i]*ker_y[j];

        if ((nvecs > 0) && (mode))
          for(idx=0;idx<kern*kern;idx++)
            vecs[nvecs][idx] -= vecs[0][idx];

        nvecs++;
      }
    }
  }

  if (nvecs != ker->nvecs)
  {
    fprintf(stderr, "make_vectors: warning: wrong number of ");
    fprintf(stderr, "kernel vectors - shouldn't happen!\n");
    fflush(stderr);
  }

  free(ker_x);
  free(ker_y);

  return;
}
