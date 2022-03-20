/*========================================================*/
/*                                                        */
/*  make_vectors.c      version 1.3.0   2005.03.03        */
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

#include <math.h>
#include <malloc.h>

#include "defs.h"
#include "funcs.h"
#include "errmess.h"

void make_vectors(double *imref, double **tab00, double **vecs, double **wxy,
                  PARAMS par)
{
        int     nvecs, nspat, nbkg, nwxy, ntot, mode,
                i, j, m, n, k, kern, idx, hnx, hny;
        double  *ker_x, *ker_y, *tmp,
                px, py, ax, ay, dx, dy;

  kern  = 2*par.kerhw + 1;

  if (!(ker_x = (double *)malloc(kern*sizeof(double)))) errmess("malloc(ker_x)");
  if (!(ker_y = (double *)malloc(kern*sizeof(double)))) errmess("malloc(ker_y)");
  if (!(tmp = (double *)malloc(par.nx*par.ny*sizeof(double))))
    errmess("malloc(tmp)");

  hnx = par.nx/2;
  hny = par.ny/2;

  dx = 1.0/(hnx - par.kerhw);
  dy = 1.0/(hny - par.kerhw);

/* vectors for background and spatial kernel variability together ! */

  for (i=par.kerhw; i<par.nx-par.kerhw; i++)
  {
    px = (i - hnx)*dx;

    for (j=par.kerhw; j<par.ny-par.kerhw; j++)
    {
      py = (j - hny)*dy;
      idx = i + par.nx*j;

      nspat = 0;
      for (m=0, ax=1.0; m<=par.sdeg; m++, ax*=px)
        for(n=0, ay=1.0; n<=par.sdeg-m; n++, ay*=py)
          tab00[nspat++][idx] = ax*ay;

      if (nspat != par.nspat)
      {
        fprintf(stderr, "make_vectors: warning: wrong number of ");
        fprintf(stderr, "spatial vectors - shouldn't happen !\n");
        fflush(stderr);
      }
    }
  }

/* now simply make pointers to appropriate elements of tab00 */

  nspat = nbkg = nwxy = 0;

  for (m=0; m<=par.sdeg; m++)
  {
    for (n=0; n<=par.sdeg-m; n++)
    {
      if ((m <= par.bdeg) && (n <= par.bdeg-m)) vecs[nbkg++] = tab00[nspat];
      if ((m <= par.wdeg) && (n <= par.wdeg-m)) wxy[nwxy++]  = tab00[nspat];
      nspat++;
    }
  }

/* the same trick for constant part of the kernel vectors */

  nvecs = nbkg;
  ntot  = nspat;

  for (k=0; k<par.ncomp; k++)
  {
    for (m=0; m<=par.deg[k]; m++)
    {
      for (n=0; n<=par.deg[k]-m; n++)
      {
        vecs[nvecs] = tab00[ntot++];

        mode = ((m/2)*2-m+1)*((n/2)*2-n+1);

        base_func(par.sig[k], m, n, ker_x, ker_y, par.kerhw, mode);
        xy_convolve(imref, vecs[nvecs], tmp, par.nx, par.ny, ker_x, ker_y,
                    par.kerhw);

/* the first vector is normalized to 1.0                     */
/* the rest must have zero sums to enforce flux conservation */

        if ((nvecs > nbkg) && (mode))
          for(idx=0; idx<par.nx*par.ny; idx++)
            vecs[nvecs][idx] -= vecs[nbkg][idx];

        nvecs++;
      }
    }
  }

  if (nvecs != par.nvecs)
  {
    fprintf(stderr, "make_vectors: warning: wrong number of ");
    fprintf(stderr, "kernel vectors - shouldn't happen !\n");
    fflush(stderr);
  }

  free(ker_x);
  free(ker_y);
  free(tmp);

  return;
}
