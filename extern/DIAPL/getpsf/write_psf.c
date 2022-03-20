/*========================================================*/
/*                                                        */
/*  write_psf.c         version 1.3.1   2005.01.28        */
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

#include "defs.h"
#include "errmess.h"

void write_psf(int vopt, char *fname, PSF_STRUCT *psf)
{
        int   i,
              sdeg, ldeg, ntot;
        FILE  *outf;

  sdeg = psf->ndeg_spat;
  ldeg = psf->ndeg_local;

  ntot = psf->ngauss*(ldeg+1)*(ldeg+2)*(sdeg+1)*(sdeg+2)/4;

  if (!(outf = fopen(fname, "wb"))) errmess(fname);

  fwrite(&psf->hw        , sizeof(int), 1, outf);
  fwrite(&psf->ndeg_spat , sizeof(int), 1, outf);
  fwrite(&psf->ndeg_local, sizeof(int), 1, outf);
  fwrite(&psf->ngauss    , sizeof(int), 1, outf);
  fwrite(&psf->recenter  , sizeof(int), 1, outf);

  fwrite(&psf->cos, sizeof(double), 1, outf);
  fwrite(&psf->sin, sizeof(double), 1, outf);
  fwrite(&psf->ax , sizeof(double), 1, outf);
  fwrite(&psf->ay , sizeof(double), 1, outf);

  fwrite(&psf->sigma_inc   , sizeof(double), 1, outf);
  fwrite(&psf->sigma_mscale, sizeof(double), 1, outf);

  fwrite(&psf->fitrad, sizeof(double), 1, outf);
  fwrite(&psf->x_orig, sizeof(double), 1, outf);
  fwrite(&psf->y_orig, sizeof(double), 1, outf);

  fwrite(psf->vec, sizeof(double), ntot, outf);

  fclose(outf);

  if (vopt > 2)
  {
    printf("PSF:\n");
    printf("====\n");
    printf("psf->hw = %d\n", psf->hw);
    printf("psf->ndeg_spat = %d\n", psf->ndeg_spat);
    printf("psf->ndeg_local = %d\n", psf->ndeg_local);
    printf("psf->ngauss = %d\n", psf->ngauss);
    printf("psf->recenter = %d\n", psf->recenter);

    printf("psf->cos = %g\n", psf->cos);
    printf("psf->sin = %g\n", psf->sin);
    printf("psf->ax = %g\n", psf->ax);
    printf("psf->ay = %g\n", psf->ay);

    printf("psf->sigma_inc = %g\n", psf->sigma_inc);
    printf("psf->sigma_mscale = %g\n", psf->sigma_mscale);

    printf("psf->fitrad = %g\n", psf->fitrad);
    printf("psf->x_orig = %g\n", psf->x_orig);
    printf("psf->y_orig = %g\n", psf->y_orig);

    for (i=0; i<ntot; i++)
      printf("  psf->vec[%d] = %g\n", i, psf->vec[i]);
    printf("====\n");
  }

  return;
}
