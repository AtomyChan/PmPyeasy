/*========================================================*/
/*                                                        */
/*  psf_stars.c         version 1.3.1   2005.01.28        */
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
#include <math.h>
#include <stdlib.h>

#include "defs.h"
#include "funcs.h"
#include "errmess.h"

#define BIG_SIGMA 100.0

void psf_stars(IMAGE *im, OBJ_LIST *objects, OBJ_LIST *psf_objects, FITPAR par)
{
        int     k, npsf, nobj, ipsf, jpsf, iobj, nbox_x, nbox_y, size_x, size_y,
                npsf_max, *nbox, xc, yc, xbox, ybox, idx, nbox_max, psfhw,
                *psf_index, *obj_index, nx, ny;
        double  flux, *mean_rat, *sigma_rat, dif, rat;
        OBJECT  *obj, *psf_obj;

  nx = im->nx;
  ny = im->ny;

  npsf = par.npsf_max;
  nobj = objects->nobj;

  nbox_x = par.nbox_x;
  nbox_y = par.nbox_y;
  psfhw  = par.psfhw;

  npsf_max = par.npsf_max;
  nbox_max = (int)ceil((double)npsf_max/(double)(nbox_x*nbox_y));

  size_x = (int)ceil((double)nx/(double)nbox_x);
  size_y = (int)ceil((double)ny/(double)nbox_y);

  if (!(mean_rat  = (double *)malloc(nbox_x*nbox_y*sizeof(double))))
    errmess("malloc(mean_rat)");
  if (!(sigma_rat = (double *)malloc(nbox_x*nbox_y*sizeof(double))))
    errmess("malloc(sigma_rat)");
  if (!(nbox = (int *)malloc(nbox_x*nbox_y*sizeof(int))))
    errmess("malloc(box)");

  for (k=0; k<nbox_x*nbox_y; k++)
  {
    mean_rat[k] = sigma_rat[k] = 0.0;
    nbox[k] = 0;
  }

  obj_index = objects->index;
  if (!(psf_index = (int *)malloc(nobj*sizeof(int))))
    errmess("malloc(psf_index)");

/* take brightest stars and equalize roughly the density over the frame */

  for (ipsf=0, iobj=nobj-1; ipsf<npsf && iobj>=0; iobj--)
  {
    obj = &objects->list[obj_index[iobj]];

    xc = obj->x;
    yc = obj->y;

    xbox = floor((double)xc/(double)size_x);
    ybox = floor((double)yc/(double)size_y);

    idx = xbox + nbox_x*ybox;

    if (nbox[idx] <= nbox_max)
    {
      if ((xc >= psfhw) && (xc < nx-psfhw) && (yc >= psfhw) && (yc < ny-psfhw))
      {
        center(im, obj);
        refine_max(im, obj);

        dif  = obj->max - obj->bg;
        flux = obj->flux;
        rat  = dif/(flux+0.01);

        if ((rat > 0.0) && (rat <= par.rat_thresh))
        {
          if ((dif > par.max_thresh) && (flux > par.min_flux))
          {
            mean_rat[idx]  += rat;
            sigma_rat[idx] += rat*rat;
            nbox[idx]++;
            obj->rat = rat;
            psf_index[ipsf++] = obj_index[iobj];
            obj->box_index = idx;
          }
        }
      }
    }
  }

  for (k=0; k<nbox_x*nbox_y; k++)
  {
    if (nbox[k] > 0)
    {
      mean_rat[k] /= nbox[k];
      if (nbox[k] >= par.min_nbox)
        sigma_rat[k] = sqrt(sigma_rat[k]/nbox[k] - mean_rat[k]*mean_rat[k]);
      else
        sigma_rat[k] = BIG_SIGMA;
    }
    else
    {
      sigma_rat[k] = mean_rat[k] = BIG_SIGMA;
    }
  }

  npsf = ipsf;

  if (!(psf_obj = (OBJECT *)malloc(npsf*sizeof(OBJECT))))
    errmess("malloc(psf_obj)");

  for (ipsf=0, jpsf=0; ipsf<npsf; ipsf++)
  {
    obj = &objects->list[psf_index[ipsf]];

    rat = obj->rat;
    idx = obj->box_index;

    xc = obj->x;
    yc = obj->y;

    if (fabs(rat-mean_rat[idx]) < par.nsig_rat*sigma_rat[idx])
    {
      if (is_isolated(im, obj, par))
      {
        copy_obj(&psf_obj[jpsf], obj);
        psf_index[jpsf] = psf_index[ipsf];
        jpsf++;
      }
    }
  }

  if (par.verbose)  printf("Found %d psf stars\n", jpsf);

  psf_objects->index = psf_index;
  psf_objects->nobj  = jpsf;
  psf_objects->list  = psf_obj;

  return;
}
