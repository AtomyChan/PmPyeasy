/*========================================================*/
/*                                                        */
/*  local_max.c         version 1.3.2   2005.04.18        */
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

void local_max(IMAGE *im, FITPAR par, OBJ_LIST *objects)
{
        int         nobj, is_max, i, j, k, l, nx, ny,
                    *obj_index;
        float       *tmp_flux;
        double      nsig2, nsig2_detect, bg, tmp_max;
        OBJECT      *obj;
        BKG_STRUCT  bkg;

  nx = im->nx;
  ny = im->ny;

  bkg.frac = par.bkg_frac;

  nobj = 0;
  if (!(obj=(OBJECT *)calloc(1,sizeof(OBJECT)))) errmess("calloc(obj)");

/* prepare some background stuff that is not changing to speed things up */
  init_bkg(im, &bkg, par.anrad1, par.anrad2);

/* loop over all pixels */

  nsig2_detect = par.nsig_detect*par.nsig_detect;

  for (i=par.psfhw; i<nx-par.psfhw; i++)
  {
    for (j=par.psfhw; j<ny-par.psfhw; j++)
    {
      tmp_max = im->pix[i+nx*j];

      is_max = 1;

      for (k=-par.maxhw; k<=par.maxhw && is_max; k++)
        for (l=-par.maxhw; l<=par.maxhw && is_max; l++)
          if (im->pix[i+k+nx*(j+l)] > tmp_max)  is_max = 0;

      if (is_max && is_peak(im, i, j, par.peakhw, par.contrast))
      {
        bg = getbkg(im, &bkg, i, j);

        nsig2 = (tmp_max - bg)*(tmp_max - bg)*im->gain/tmp_max;

        if (!(obj = (OBJECT *)realloc(obj, (nobj+1)*sizeof(OBJECT))))
          errmess("realloc(obj)");

        if ((nsig2 >= nsig2_detect) && (bg > 0.0))
        {
          obj[nobj].x = i;
          obj[nobj].y = j;
          obj[nobj].max = tmp_max;
          obj[nobj].bg = bg;
          obj[nobj].flux = silly_phot(im, &obj[nobj], (int)(par.aprad+0.5));
          obj[nobj].norm = obj[nobj].flux;

          nobj++;
        }
      }
/* end of loop over all pixels */
    }
  }

  if (par.verbose)  printf("Found %d good local maxima\n", nobj);

  objects->list = obj;
  objects->nobj = nobj;

  if (!(tmp_flux=(float *)malloc(nobj*sizeof(float))))
    errmess("malloc(tmp_flux)");

  if (!(obj_index = (int *)malloc(nobj*sizeof(int))))
    errmess("malloc(obj_index)");

  for (i=0; i<nobj; i++)  tmp_flux[i] = obj[i].flux;

  quick_sort(tmp_flux, obj_index, nobj);

  objects->index = obj_index;

  free(tmp_flux);

  return;
}
