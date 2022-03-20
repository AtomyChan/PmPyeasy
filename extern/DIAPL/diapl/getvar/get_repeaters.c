/*========================================================*/
/*                                                        */
/*  get_repeaters.c     version 1.2.3   2005.01.24        */
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
#include <stdlib.h>

#include "defs.h"
#include "funcs.h"
#include "errmess.h"

#define BIG_double 1.0e32

void get_repeaters(double **data, double **noise, double *varim1, double *varim2,
                   PAR_STRUCT *par, int nim)
{
        int   i, j, k, l, nmarg, nx, ny, *index, iim, ngood, ncons_var1,
              npts_var2, idx1, idx2, hw, npix,
              nhigh1, nhigh2, ncons1, ncons2, npts1, npts2, ncmax1, ncmax2;
        double gain, sat_level, base_flux, *tab, *sig,
              nsig, nsig_var1, nsig_var2,
              chi2_1, chi2_2, ratio, lim_ratio, tmp, tmp_sig,
              pix, pix_sig;

  nx = par->nx0;
  ny = par->ny0;
  hw = par->filter_hw;

  nmarg = (int)par->bad_margin - hw - 3;

  if (nmarg < hw) nmarg = hw;

  gain  = par->gain;

  nsig_var1  = par->nsig_var1;
  nsig_var2  = par->nsig_var2;
  ncons_var1 = par->ncons_var1;
  npts_var2  = par->npts_var2;
  lim_ratio  = par->lim_ratio;

  sat_level  = par->sat_level;

  if (!(tab   = (double *)malloc(nim*sizeof(double)))) errmess("malloc(tab)");
  if (!(sig   = (double *)malloc(nim*sizeof(double)))) errmess("malloc(sig)");
  if (!(index = (int   *)malloc(nim*sizeof(int)))) errmess("malloc(index)");

  for (i=nmarg; i<nx-nmarg; i++)
  {
    for (j=nmarg; j<ny-nmarg; j++)
    {
      idx1 = i + nx*j;

      for (iim=0; iim<nim; iim++)
      {
        pix  = pix_sig = 0.0;
        npix = 0;

        for (k=-hw; k<=hw; k++)
        {
          for (l=-hw; l<=hw; l++)
          {
            idx2 = idx1 + k + nx*l;
            tmp     =  data[iim][idx2];
            tmp_sig = noise[iim][idx2];

            if ((tmp < sat_level) && (tmp_sig < sat_level))
            {
              pix += tmp;
              pix_sig += tmp_sig;
              npix++;
            }
          }
        }

        if (npix > 0)
        {
          tab[iim] = pix/npix;
          sig[iim] = pix_sig/(npix*npix);
        }
        else
        {
          tab[iim] = BIG_double;
          sig[iim] = BIG_double;
        }

        if ((sig[iim] > 0.0) && (sig[iim] < BIG_double))
          sig[iim] = sqrt(sig[iim]/gain);
        else
          sig[iim] = BIG_double;
      }

      quick_sort(tab, index, nim);

      ngood = 0;

      while ((ngood < nim) && (tab[index[ngood]] < BIG_double))  ngood++;

      base_flux = tab[index[ngood/2]];

      ncons1 = ncons2 = 0;
      nhigh1 = nhigh2 = 0;
      npts1  = npts2  = 0;
      ncmax1 = ncmax2 = 0;

      chi2_1 = chi2_2 = 0.0;

      if (ngood >= nim/2)
      {
        for (iim=0; iim<nim; iim++)
        {
          if ((tab[iim] < BIG_double) && (sig[iim] < BIG_double))
          {
            nsig  = (tab[iim] - base_flux)/sig[iim];

            if (nsig >= nsig_var1)
            {
              ncons1++;
              npts1++;
              chi2_1 += tab[iim] - base_flux;
            }
            else  ncons1=0;

            if (nsig <= -nsig_var1)
            {
              ncons2++;
              npts2++;
              chi2_2 += tab[iim] - base_flux;
            }
            else  ncons2=0;

            if (nsig >= nsig_var2)  nhigh1++;
            if (nsig <= -nsig_var2) nhigh2++;

            if (ncons1 > ncmax1)  ncmax1 = ncons1;
            if (ncons2 > ncmax2)  ncmax2 = ncons2;
          }
        }
      }

      varim1[idx1] = varim2[idx1] = 0.0;

      if (npts2 > 0)  ratio = (double)npts1/(double)npts2;
      else            ratio = BIG_double;

      if ((nhigh1 >= npts_var2) || (nhigh2 >= npts_var2) ||
          (ncmax1 >= ncons_var1) || (ncmax2 >= ncons_var1))
      {
        if ((ratio < lim_ratio) && (ratio > 1.0/lim_ratio))
          varim1[idx1] = (fabs(chi2_1) + fabs(chi2_2))/(npts1+npts2);
        else if (ratio > lim_ratio)
          varim2[idx1] = fabs(chi2_1)/npts1;
        else
          varim2[idx1] = fabs(chi2_2)/npts2;
      }
    }
  }

  free(tab);
  free(sig);
  free(index);

  return;
}
