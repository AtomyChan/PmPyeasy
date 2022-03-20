/*========================================================*/
/*                                                        */
/*  read_psf.c          version 1.3.0   2005.04.20        */
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
#include "errmess.h"

void read_psf(char *fname, PSF_STRUCT *psf, int verbose)
{
        int   i,
              sdeg, ldeg, ntot, nloc;
        FILE  *inpf;

  if (!(inpf = fopen(fname, "r"))) errmess(fname);

  fread(&psf->hw        , sizeof(int), 1, inpf);
  fread(&psf->ndeg_spat , sizeof(int), 1, inpf);
  fread(&psf->ndeg_local, sizeof(int), 1, inpf);
  fread(&psf->ngauss    , sizeof(int), 1, inpf);
  fread(&psf->recenter  , sizeof(int), 1, inpf);

  fread(&psf->cos , sizeof(double), 1, inpf);
  fread(&psf->sin , sizeof(double), 1, inpf);
  fread(&psf->ax  , sizeof(double), 1, inpf);
  fread(&psf->ay  , sizeof(double), 1, inpf);

  fread(&psf->sigma_inc   , sizeof(double), 1, inpf);
  fread(&psf->sigma_mscale, sizeof(double), 1, inpf);

  fread(&psf->fitrad , sizeof(double), 1, inpf);
  fread(&psf->x_orig , sizeof(double), 1, inpf);
  fread(&psf->y_orig , sizeof(double), 1, inpf);

  sdeg = psf->ndeg_spat;
  ldeg = psf->ndeg_local;

  nloc = psf->ngauss*(ldeg+1)*(ldeg+2)/2;
  ntot =        nloc*(sdeg+1)*(sdeg+2)/2;

  if (!(psf->vec = (double *)malloc(ntot*sizeof(double))))
    errmess("malloc(psf->vec)");
  if (!(psf->local_vec = (double *)malloc(nloc*sizeof(double))))
    errmess("malloc(psf->local_vec)");

  fread(psf->vec , sizeof(double), ntot, inpf);

  fclose(inpf);

  if (verbose > 2)
  {
    printf("Reading: %s\n", fname);
    printf("-----------\n");
    printf("psf->hw        = %d\n", psf->hw);
    printf("psf->ndeg_spat = %d\n", psf->ndeg_spat);
    printf("psf->ndeg_local= %d\n", psf->ndeg_local);
    printf("psf->ngauss    = %d\n", psf->ngauss);
    printf("psf->recenter  = %d\n", psf->recenter);

    printf("psf->cos = %g\n", psf->cos);
    printf("psf->sin = %g\n", psf->sin);
    printf("psf->ax  = %g\n", psf->ax);
    printf("psf->ay  = %g\n", psf->ay);

    printf("psf->sigma_inc   = %g\n", psf->sigma_inc);
    printf("psf->sigma_mscale= %g\n", psf->sigma_mscale);

    printf("psf->fitrad = %g\n", psf->fitrad);
    printf("psf->x_orig = %g\n", psf->x_orig);
    printf("psf->y_orig = %g\n", psf->y_orig);

    printf("psf->ndeg_spat = %d\n", psf->ndeg_spat);
    printf("psf->ndeg_local= %d\n", psf->ndeg_local);
    printf("nloc= %d\n", nloc);
    printf("ntot= %d\n", ntot);

    if (verbose > 3)
      for (i=0; i<ntot; i++) printf("  psf->vec[%d]= %g\n", i, psf->vec[i]);
  }

  return;
}
