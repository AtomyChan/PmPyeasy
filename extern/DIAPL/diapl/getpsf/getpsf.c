/*========================================================*/
/*                                                        */
/*  getpsf.c            version 1.4.0   2005.04.18        */
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
/* Program determines spatially variable PSF for an image */
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
#include "pfitsin1.h"

int main(int argc, char *argv[])
{
        char        **header, *imname, *parfname,
                    *instrname,
                    *outfname;
        int         ntot, ncomp, k, nstars, headlen, i;
        FITPAR      par;
        OBJ_LIST    objects, psf_objects;
        PSF_STRUCT  psf;
        IMAGE       im;

/* IO stuff */
  if (argc != 5)
  {
    printf("\n\tUSAGE: getpsf  parameter_file instrument_file ");
    printf("input_image output_file\n");
    exit(1);
  }

  parfname = argv[1];
  instrname= argv[2];
  imname   = argv[3];
  outfname = argv[4];

/* get parameters */
  get_params(parfname, instrname, &im, &par, &psf);

  if (par.verbose > 2)
  {
    printf("imname=   %s\n", imname);
    printf("outfname= %s\n", outfname);
  }

  psf.fitrad = par.fitrad;
  psf.ngauss = par.ngauss;
  psf.hw     = par.psfhw;
  psf.gain   = im.gain;
  psf.x_orig = im.nx/2;
  psf.y_orig = im.ny/2;

/* read image */
  im.pix=read_FITS_2D1file(imname, 's', &headlen, &header, &im.nx, &im.ny);
  for (i=0; i<headlen; i++) free(header[i]);
  free(header);

  ncomp = par.ngauss*(par.ndeg_local+1)*(par.ndeg_local+2)/2;
  ntot  = ncomp*(par.ndeg_spat+1)*(par.ndeg_spat+2)/2;

  if (!(psf.vec=(double *)malloc(ntot*sizeof(double))))
    errmess("malloc(psf.vec)");
  if (!(psf.local_vec=(double *)malloc(ncomp*sizeof(double))))
    errmess("malloc(psf.local)");

/* find good stars for PSF fit and get relevant pixels */
  if (par.verbose)  printf("\n");

  local_max(&im, par, &objects);

  psf_stars(&im, &objects, &psf_objects, par);

  psf_bitmaps(&im, &psf_objects, &psf);

/* initial fit */
  psf.ndeg_spat  = 0;
  psf.ndeg_local = par.ndeg_local;

  if (par.verbose)
  {
    printf("\nInitial fit (constant PSF):\n\n");
    printf("iter\t      ax\t      ay\t    cos(fi)\t    sin(fi)  \n");
  }

  for (k=1; k<=par.niter_init; k++)
  {
    fitpsf(&psf_objects, &psf);

    get_gaussian(&psf);

    local_fit(&psf_objects, &psf);

    if (par.verbose)
      printf("%2d\t%10.3f\t%10.3f\t%10.3f\t%10.3f\n",
              k, psf.ax, psf.ay, psf.cos, psf.sin);
  }

/* final fit with clipping */
  psf.ndeg_spat  = par.ndeg_spat;
  psf.ndeg_local = par.ndeg_local;

  if (par.verbose)
  {
    printf("\nFinal fit (variable PSF):\n\n");
    printf("iter\t     N(good)    N(clip) \n");
  }

  for (k=1; k<=par.niter; k++)
  {
    clip_stars(&psf_objects, par.nsig_clip, &nstars);

    fitpsf(&psf_objects, &psf);

    local_fit(&psf_objects, &psf);

    if (par.verbose)
      printf("%2d\t      %3d\t %3d\n", k, nstars, psf_objects.nobj-nstars);
  }

/* output to a binary file */
  write_psf(par.verbose, outfname, &psf);

  if (par.verbose)
    printf("\n\nPSF coefficients written to %s\n", outfname);

  return(0);
}
