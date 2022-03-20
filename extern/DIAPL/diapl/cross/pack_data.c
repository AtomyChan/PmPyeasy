/*========================================================*/
/*                                                        */
/*  pack_data.c         version 1.3.5   2005.04.26        */
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

#include <string.h>

#include "defs.h"
#include "errmess.h"

void pack_data(double *im, double *data, PARAMS par)
{
        long  *hist, idx,
              hist_size,
              npix;
        int   i, j, k, l;
        double bkg, sum, frac, norm, pix;

  norm=par.x_nbin*par.y_nbin;

  for (j=0; j<par.ny; j++)
  {
    for (i=0; i<par.nx; i++)
    {
      pix=0.0;

      for (l=0; l<par.y_nbin; l++)
      {
	for (k=0; k<par.x_nbin; k++)
        {
	  idx = i*par.x_nbin + k + par.nx0*(j*par.y_nbin + l);

	  if      (im[idx] > par.sat_level) pix+=par.sat_level;
          else if (im[idx] < par.min_level) pix+=par.min_level;
	  else                              pix+=im[idx];
	}
      }

      idx=2*(i + par.nx*j);

      data[idx]=pix/norm;
      data[idx+1]=0.0;
    }
  }

  hist_size=(long)(par.sat_level-par.min_level+1);
  if (!(hist=(long *)calloc(hist_size, sizeof(long))))
    errmess("calloc(hist)");
  memset(hist, 0, hist_size*sizeof(long));

  npix=0L;
  for(i=0; i<2*par.nx*par.ny; i+=2)
  {
    idx=(long)(data[i]-par.min_level);
    if ((idx >= 0) && (idx < hist_size))
    {
      hist[idx]++;
      npix++;
    }
    else
    {
      printf("ERROR! pack_data(): idx= %ld  hist_size= %ld\n", idx, hist_size);
      exit(22);
    }
  }

  sum=par.bkg_frac*npix;

  frac=0.0;
  idx=0L;
  while (frac < sum) frac+=hist[idx++];

  free(hist);

  bkg=idx+par.min_level-0.5;
  if (par.verbose > 1) printf("background= %g\n", bkg);

  for(i=0; i<2*par.nx*par.ny; i+=2) data[i]-=bkg;

  return;
}
