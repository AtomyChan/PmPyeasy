/*========================================================*/
/*                                                        */
/*  get_params.c        version 1.5.1   2006.03.31        */
/*  (getpsf)                                              */
/*                                                        */
/*  Copyright (C) 2006 by Wojtek Pych, CAMK PAN           */
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
#include <string.h>

#include "defs.h"
#include "errmess.h"

void get_params(char *parfname, char *instrname, IMAGE *im, FITPAR *par,
                PSF_STRUCT *psf)
{
        char  line[201],
              key[200],
              val[200];
        int   i;
        FILE *inpf;

  if (!(inpf = fopen(parfname, "r"))) errmess(parfname);

  for (i=0; !feof(inpf); i++)
  {
    fgets(line, 200, inpf);
    if (feof(inpf)) break;

    sscanf(line, "%s = %s", key, val);

    if (!strcasecmp(key, "END")) break;
    else if (!strcasecmp(key, "NBOX_X"))     sscanf(val, "%d", &par->nbox_x);
    else if (!strcasecmp(key, "NBOX_Y"))     sscanf(val, "%d", &par->nbox_y);
    else if (!strcasecmp(key, "NDEG_SPAT"))  sscanf(val, "%d", &par->ndeg_spat);
    else if (!strcasecmp(key, "NDEG_LOCAL")) sscanf(val, "%d", &par->ndeg_local);
    else if (!strcasecmp(key, "NGAUSS"))     sscanf(val, "%d", &par->ngauss);
    else if (!strcasecmp(key, "NPSF_MAX"))   sscanf(val, "%d", &par->npsf_max);
    else if (!strcasecmp(key, "MIN_NBOX"))   sscanf(val, "%d", &par->min_nbox);
    else if (!strcasecmp(key, "MIN_FLUX"))   sscanf(val, "%lf", &par->min_flux);
    else if (!strcasecmp(key, "MAX_THRESH")) sscanf(val, "%lf", &par->max_thresh);
    else if (!strcasecmp(key, "RAT_THRESH")) sscanf(val, "%lf", &par->rat_thresh);
    else if (!strcasecmp(key, "NSIG_DETECT")) sscanf(val, "%lf", &par->nsig_detect);
    else if (!strcasecmp(key, "CONTRAST"))   sscanf(val, "%lf", &par->contrast);
    else if (!strcasecmp(key, "NSIG_RAT"))   sscanf(val, "%lf", &par->nsig_rat);
    else if (!strcasecmp(key, "PSFHW"))      sscanf(val, "%d", &par->psfhw);
    else if (!strcasecmp(key, "MAXHW"))      sscanf(val, "%d", &par->maxhw);
    else if (!strcasecmp(key, "PEAKHW"))     sscanf(val, "%d", &par->peakhw);
    else if (!strcasecmp(key, "ISOHW"))      sscanf(val, "%d", &par->isohw);
    else if (!strcasecmp(key, "ISO_OFF"))    sscanf(val, "%lf", &par->iso_off);
    else if (!strcasecmp(key, "ISO_SLO"))    sscanf(val, "%lf", &par->iso_slo);
    else if (!strcasecmp(key, "FITRAD"))     sscanf(val, "%lf", &par->fitrad);
    else if (!strcasecmp(key, "APRAD"))      sscanf(val, "%lf", &par->aprad);
    else if (!strcasecmp(key, "ANRAD1"))     sscanf(val, "%lf", &par->anrad1);
    else if (!strcasecmp(key, "ANRAD2"))     sscanf(val, "%lf", &par->anrad2);
    else if (!strcasecmp(key, "BKG_FRAC"))   sscanf(val, "%lf", &par->bkg_frac);
    else if (!strcasecmp(key, "NSIG_CLIP"))  sscanf(val, "%lf", &par->nsig_clip);
    else if (!strcasecmp(key, "NITER_INIT")) sscanf(val, "%d", &par->niter_init);
    else if (!strcasecmp(key, "NITER"))      sscanf(val, "%d", &par->niter);
    else if (!strcasecmp(key, "RECENTER"))   sscanf(val, "%d", &psf->recenter);
    else if (!strcasecmp(key, "PSF_COS"))    sscanf(val, "%lf", &psf->cos);
    else if (!strcasecmp(key, "PSF_SIN"))    sscanf(val, "%lf", &psf->sin);
    else if (!strcasecmp(key, "PSF_AX"))     sscanf(val, "%lf", &psf->ax);
    else if (!strcasecmp(key, "PSF_AY"))     sscanf(val, "%lf", &psf->ay);
    else if (!strcasecmp(key, "SIGMA_INC"))  sscanf(val, "%lf", &psf->sigma_inc);
    else if (!strcasecmp(key, "SIGMA_MSCALE")) sscanf(val, "%lf", &psf->sigma_mscale);
    else if (!strcasecmp(key, "VERBOSE"))    sscanf(val, "%d", &par->verbose);
    else printf("WARNING: Unknown parameter %s\n", key);
  }

  fclose(inpf);

  if (par->verbose > 1)
  {
    printf("Reading %s\n", parfname);
    printf("----------\n");
    printf("NBOX_X        = %d\n", par->nbox_x);
    printf("NBOX_Y        = %d\n", par->nbox_y);
    printf("NDEG_SPAT     = %d\n", par->ndeg_spat);
    printf("NDEG_LOCAL    = %d\n", par->ndeg_local);
    printf("NGAUSS        = %d\n", par->ngauss);
    printf("NPSF_MAX      = %d\n", par->npsf_max);
    printf("MIN_NBOX      = %d\n", par->min_nbox);
    printf("MIN_FLUX      = %g\n", par->min_flux);
    printf("MAX_THRESH    = %g\n", par->max_thresh);
    printf("RAT_THRESH    = %g\n", par->rat_thresh);
    printf("NSIG_DETECT   = %g\n", par->nsig_detect);
    printf("CONTRAST      = %g\n", par->contrast);
    printf("NSIG_RAT      = %g\n", par->nsig_rat);
    printf("PSFHW         = %d\n", par->psfhw);
    printf("MAXHW         = %d\n", par->maxhw);
    printf("PEAKHW        = %d\n", par->peakhw);
    printf("ISOHW         = %d\n", par->isohw);
    printf("ISO_OFF       = %g\n", par->iso_off);
    printf("ISO_SLO       = %g\n", par->iso_slo);
    printf("FITRAD        = %g\n", par->fitrad);
    printf("APRAD         = %g\n", par->aprad);
    printf("ANRAD1        = %g\n", par->anrad1);
    printf("ANRAD2        = %g\n", par->anrad2);
    printf("BKG_FRAC      = %g\n", par->bkg_frac);
    printf("NSIG_CLIP     = %g\n", par->nsig_clip);
    printf("NITER_INIT    = %d\n", par->niter_init);
    printf("NITER         = %d\n", par->niter);
    printf("RECENTER      = %d\n", psf->recenter);
    printf("PSF_COS       = %g\n", psf->cos);
    printf("PSF_SIN       = %g\n", psf->sin);
    printf("PSF_AX        = %g\n", psf->ax);
    printf("PSF_AY        = %g\n", psf->ay);
    printf("SIGMA_INC     = %g\n", psf->sigma_inc);
    printf("SIGMA_MSCALE  = %g\n", psf->sigma_mscale);
    printf("VERBOSE       = %d\n", par->verbose);
    printf("----------\n");
  }


  if (!(inpf = fopen(instrname, "r"))) errmess(instrname);

  for (i=0; !feof(inpf); i++)
  {
    fgets(line, 200, inpf);
    if (feof(inpf)) break;

    sscanf(line, "%s = %s", key, val);

    if (!strcasecmp(key, "END")) break;
    else if (!strcasecmp(key, "GAIN"))      sscanf(val, "%lf", &im->gain);
    else if (!strcasecmp(key, "SAT_LEVEL")) sscanf(val, "%lf", &im->sat_level);
    else if (!strcasecmp(key, "MIN_LEVEL")) sscanf(val, "%lf", &im->min_level);
    else if (!strcasecmp(key, "NX"))        continue;
    else if (!strcasecmp(key, "NY"))        continue;
    else if (!strcasecmp(key, "MARG"))      continue;
    else if (!strcasecmp(key, "FITS"))      continue;
    else printf("WARNING: Unknown parameter %s\n", key);
  }

  fclose(inpf);

  if (par->verbose > 1)
  {
    printf("Reading %s\n", instrname);
    printf("----------\n");
    printf("GAIN          = %g\n", im->gain);
    printf("SAT_LEVEL     = %g\n", im->sat_level);
    printf("MIN_LEVEL     = %g\n", im->min_level);
    printf("----------\n");
  }

  return;
}
/*** END ***/
