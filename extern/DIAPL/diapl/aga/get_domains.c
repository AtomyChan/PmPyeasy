/*========================================================*/
/*                                                        */
/*  get_domains.c       version 1.3.0   2005.03.03        */
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

void get_domains(double *im, unsigned short *mask0, DOM_TABS *domp,
                 PARAMS par, int *ndom)
{
        int   i, j, nx_sector, ny_sector, xs, ys,
              i_max=0, j_max=0, idx, idom, idom_x, idom_y, nmarg;
        double pixmax, pix;

  nmarg = par.kerhw + par.domhw;

  idom = 0;

  if (par.domain_mode)
  {
/* domains around bright stars */
    nx_sector = (int)ceil((par.nx-2*nmarg)/(double)par.ndom_x);
    ny_sector = (int)ceil((par.ny-2*nmarg)/(double)par.ndom_y);

    for (idom_x=0; idom_x<par.ndom_x; idom_x++)
    {
      xs = nmarg + idom_x*nx_sector;

      for (idom_y=0; idom_y<par.ndom_y; idom_y++)
      {
        ys = nmarg + idom_y*ny_sector;

        pixmax = 0.0;

        for (i=xs; i<xs+nx_sector && i<par.nx-nmarg; i++)
        {
          for (j=ys; j<ys+ny_sector && j<par.ny-nmarg; j++)
          {
            idx = i + par.nx*j;

            if (max(i, j, par.mohw, im, par.nx) && mask0[idx])
            {
              pix = im[idx];
              if (pix > pixmax)
              {
                pixmax = pix;
                i_max  = i;
                j_max  = j;
              }
            }
          }
        }

        if (pixmax > par.dom_thresh)
        {
          domp[idom].xc = i_max;
          domp[idom].yc = j_max;
          domp[idom].reject = RJCT_NONE;
          domp[idom].chi2_n = CHI2_INIT;
          idom++;
        }
      }
    }

    *ndom = idom;
  }
  else
  {
/* domains uniformly distributed in the image */

    nx_sector = (par.nx-2*nmarg-1)/(par.ndom_x-1);
    ny_sector = (par.ny-2*nmarg-1)/(par.ndom_y-1);

    for (idom_x=0; idom_x<par.ndom_x; idom_x++)
    {
      xs = nmarg + idom_x*nx_sector;

      for (idom_y=0; idom_y<par.ndom_y; idom_y++)
      {
        ys = nmarg + idom_y*ny_sector;

        if ((xs < par.nx-nmarg) && (ys < par.ny-nmarg))
        {
          domp[idom].xc = xs;
          domp[idom].yc = ys;
          domp[idom].reject = RJCT_NONE;
          domp[idom].chi2_n = CHI2_INIT;
          idom++;
        }
      }
    }

    *ndom = idom;
  }

  if (par.verbose > VERB_HIGH)
  {
    printf("%d domains:\n", *ndom);
    for(idom=0; idom<*ndom; idom++)
      printf("%3d\tx= %3d\ty= %3d\n", idom, domp[idom].xc, domp[idom].yc);
  }

  return;
}
