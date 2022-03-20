/*========================================================*/
/*                                                        */
/*  center_stars.c      version 1.2.3   2005.01.24        */
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

void center_stars(STAR *obj, int nobj, double **data, double **noise,
                  PAR_STRUCT *par, int nim)
{
        int    i, j, iim, nx, ny, iobj, nframes, size, sat_test, nsig_bkg,
               cx, cy, hw, sx, sy;
        double  gain, nsig_var1, x, y, aprad, anrad1, anrad2, sat_level,
               bg, flux, ratio, sigma, max, pix, id_rad;
        STAR   *objp;

  nsig_var1 = par->nsig_var1;
  sat_level = par->sat_level;
  nsig_bkg  = (int)(par->nsig_bkg);

  anrad1 = par->anrad1;
  anrad2 = par->anrad2;
  aprad  = par->aprad;
  id_rad = par->id_rad;
  gain   = par->gain;

  nx   = par->nx0;
  ny   = par->ny0;

  hw   = par->center_hw;
  size = 2*hw + 1;

  for (iobj=0; iobj<nobj; iobj++)
  {
    objp = &obj[iobj];

    cx = objp->cx;
    cy = objp->cy;

    nframes = 0;

    if (!(objp->bmp = (double *)malloc(size*size*sizeof(double))))
      errmess("malloc(objp->bmp)");

    for (i=0; i<size*size; i++) objp->bmp[i] = 0.0;

    for (iim=0; iim<nim; iim++)
    {
      sat_test = 1;

      bg = 0.0;

      flux=aperphot(data[iim], nx, ny, cx, cy, aprad, bg, &sat_test,
                    sat_level);
      sigma=aperphot(noise[iim], nx, ny, cx, cy, aprad, bg, &sat_test,
                     sat_level);

      if (sat_test && (sigma > 0.0))
      {
        ratio = flux/sqrt(sigma/gain);

        if ((ratio >= nsig_var1) || (ratio <= -nsig_var1))
        {
          nframes++;

          for (i=-hw; i<=hw; i++)
            for( j=-hw; j<=hw; j++)
              objp->bmp[i+hw+size*(j+hw)] += ratio*data[iim][cx+i+nx*(cy+j)];
        }
      }
    }

    if (nframes > 0)
    {
      max = objp->bmp[hw+size*hw];

      sx = sy = hw;

      for (i=-1; i<=1; i++)
      {
        for (j=-1; j<=1; j++)
        {
          pix = objp->bmp[hw+i+size*(hw+j)];

          if (pix > max)
          {
            max = pix;
            sx  = hw + i;
            sy  = hw + j;
          }
        }
      }

      centroid(objp->bmp, size, size, sx, sy, &x, &y);

      if ((fabs(x-hw) <= id_rad) || (fabs(y-hw) <= id_rad))
      {
        objp->x = cx + x - (double)hw;
        objp->y = cy + y - (double)hw;
        objp->nframes = nframes;
      }
      else
      {
        objp->x = -1.0;
        objp->y = -1.0;
        objp->nframes = 0;
      }
    }
    else
    {
      objp->x = -1.0;
      objp->y = -1.0;
      objp->nframes = 0;
    }
  }

  return;
}
