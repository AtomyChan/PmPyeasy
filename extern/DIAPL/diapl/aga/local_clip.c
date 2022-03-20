/*========================================================*/
/*                                                        */
/*  local_clip.c        version 1.3.0   2005.03.03        */
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

#include <stdlib.h>

#include "defs.h"
#include "funcs.h"
#include "errmess.h"

int local_clip(double *im, double *imref, unsigned short *mask1, double **vecs,
               double **wxy, double *exp_vec, DOM_TABS *domp, PARAMS par,
               double *chi2_n, int *npix)
{
        int     i, j, m, n, idx, nclip, idom, xc, yc, loc_npix, loc_nclip;
        double  pixdif, pixdev, chi2, n_sig2, *loc_ker, loc_chi2, dom_area,
                **mat, *vec;

  if (!(loc_ker = (double *)malloc(par.nvecs*sizeof(double))))
    errmess("malloc(loc_ker)");

  n_sig2  = (double)(par.n_sig*par.n_sig);

  dom_area  = 2*par.domhw + 1;
  dom_area *= dom_area;

  chi2  = 0.0;
  *npix = 0;
  nclip = 0;

  for (idom=0; idom<par.ndom; idom++)
  {
    loc_chi2 = 0.0;
    loc_npix = loc_nclip = 0;

    mat = domp[idom].mat;
    vec = domp[idom].vec;
    xc  = domp[idom].xc;
    yc  = domp[idom].yc;

    domp[idom].reject = RJCT_NONE;

    local_solution(xc, yc, wxy, exp_vec, loc_ker, par);

    for (j=yc-par.domhw; j<=yc+par.domhw; j++)
    {
      for (i=xc-par.domhw; i<=xc+par.domhw; i++)
      {
        idx = i + par.nx*j;
        pixdif = (double)im[idx] - loc_ker[0];

        for (m=1; m<par.nvecs; m++) pixdif -= loc_ker[m]*vecs[m][idx];

        if (mask1[idx])
        {
          pixdev = pixdif*pixdif/im[idx];

          if (pixdev < n_sig2*(*chi2_n))
          {
            loc_chi2 += pixdev;
            (*npix)++;
            loc_npix++;
          }
          else
          {
            for (m=0; m<par.nvecs; m++)
            {
              vec[m] -= im[idx]*vecs[m][idx]*imref[idx];

              for (n=0; n<=m; n++)
                mat[m][n] -= vecs[m][idx]*vecs[n][idx]*imref[idx];
            }
            nclip++;
            loc_nclip++;
          }
        }
      }
    }

    if (loc_npix > 0)
    {
      domp[idom].chi2_n = loc_chi2/loc_npix;
      chi2 += loc_chi2;
    }
    else  domp[idom].chi2_n = CHI2_INIT;

    if ((double)loc_npix < par.min_area_dom * dom_area)
      domp[idom].reject = RJCT_HIGH;

    for (m=0; m<par.nvecs; m++)
      for (n=0; n<m; n++)
        mat[n][m] = mat[m][n];
  }

  if (npix > 0) *chi2_n = chi2/(*npix);
  else          *chi2_n = CHI2_INIT;

  free(loc_ker);

  return (nclip);
}
