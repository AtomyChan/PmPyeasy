/*========================================================*/
/*                                                        */
/*  clip_domains.c      version 1.3.0   2005.03.03        */
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

#include "defs.h"
#include "funcs.h"

int clip_domains(DOM_TABS *domp, double *sort_buf, int *indx,
                 int iter, double *chi2_n, PARAMS par)
{
        int     nkeep, idom, n_sig;
        double  median, sig, sigma;

  n_sig  = par.n_sig_dom;

  nkeep = 0;
  for (idom=0; idom<par.ndom; idom++)
  {
    if (domp[idom].reject != RJCT_NONE) continue;
    sort_buf[nkeep++] = domp[idom].chi2_n;
  }

  quick_sort(sort_buf, indx, nkeep);

  median = sort_buf[indx[nkeep/2]];

  sigma = 0.0;
  for (idom=0; idom<par.ndom; idom++)
  {
    if (domp[idom].reject != RJCT_NONE) continue;
    sig    = domp[idom].chi2_n - median;
    sigma += sig*sig;
  }

  sigma /= nkeep;

  nkeep = 0;
  for (idom=0; idom<par.ndom; idom++)
  {
    sig = domp[idom].chi2_n - median;

    if (domp[idom].reject < RJCT_HIGH)
    {
      domp[idom].reject = RJCT_NONE;

      if (sig > n_sig*sqrt(sigma))  domp[idom].reject = RJCT_LOW;
      else                          nkeep++;
    }
  }
  *chi2_n = median;

  if (par.verbose >= VERB_MED)
    printf("%6d\t%14.4f\t%14.4f\t%7d\n",
            iter, par.gain*median, par.gain*sigma, nkeep);

  if (par.verbose >= VERB_HIGH)
  {
    for (idom=0; idom<par.ndom; idom++)
    {
      if (domp[idom].reject != RJCT_NONE)
      {
        printf("rejected domain %4d:\tx=%3d, \ty=%3d,",
                                idom, domp[idom].xc, domp[idom].yc);
        printf("\tchi2_n=%6.3f,\tstatus=%1d\n",
                         domp[idom].chi2_n*par.gain, domp[idom].reject);
      }
    }
  }

  return (nkeep);
}
