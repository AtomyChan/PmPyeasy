/*========================================================*/
/*                                                        */
/*  get_params.c        version 1.9.3   2006.10.25        */
/*  (mstack)                                              */
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
#include <strings.h>

#include "defs.h"
#include "errmess.h"

void get_params(char *fname, char *instrname, PAR_STRUCT *par)
{
        char  line[201],
              key[200],
              val[200];
        int   snx, sny, smarg;
        FILE  *inpf;

  if (!(inpf = fopen(fname, "r"))) errmess(fname);

  while (!feof(inpf))
  {
    if (fgets(line, 200, inpf) == NULL) errmess(fname);
    if (feof(inpf)) break;

    sscanf(line, "%s = %s", key, val);

    if (!strcasecmp(key, "END")) break;
    else if (!strcasecmp(key, "NY"))          sscanf(val, "%d", &par->ny);
    else if (!strcasecmp(key, "HIST_NBIN"))   sscanf(val, "%d", &par->nbin);
    else if (!strcasecmp(key, "HIST_LOW"))    sscanf(val, "%f", &par->low);
    else if (!strcasecmp(key, "HIST_HIGH"))   sscanf(val, "%f", &par->high);
    else if (!strcasecmp(key, "NBIN_SMOOTH")) sscanf(val, "%d", &par->nsmooth);
    else if (!strcasecmp(key, "BAD_VALUE"))   sscanf(val, "%f", &par->bad_value);
    else if (!strcasecmp(key, "BKG_FRAC"))    sscanf(val, "%f", &par->bkg_frac);
    else if (!strcasecmp(key, "THRESHOLD"))   sscanf(val, "%f", &par->threshold);
    else if (!strcasecmp(key, "MIN_SCALE"))   sscanf(val, "%f", &par->min_scale);
    else if (!strcasecmp(key, "MAX_SCALE"))   sscanf(val, "%f", &par->max_scale);
    else if (!strcasecmp(key, "MIN_NGOOD"))   sscanf(val, "%d", &par->min_ngood);
    else if (!strcasecmp(key, "N_REJ_LOW"))   sscanf(val, "%d", &par->nlow);
    else if (!strcasecmp(key, "N_REJ_HIGH"))  sscanf(val, "%d", &par->nhigh);
    else if (!strcasecmp(key, "MEDIAN"))      sscanf(val, "%d", &par->median);
    else if (!strcasecmp(key, "VERBOSE"))     sscanf(val, "%d", &par->verbose);
    else printf("WARNING: Unknown parameter %s\n", key);
  }

  fclose(inpf);

  if (par->verbose > 1)
  {
    printf("Reading %s:\n", fname);
    printf("-----------\n");
    printf("NY          = %d\n", par->ny);
    printf("HIST_NBIN   = %d\n", par->nbin);
    printf("HIST_LOW    = %g\n", par->low);
    printf("HIST_HIGH   = %g\n", par->high);
    printf("NBIN_SMOOTH = %d\n", par->nsmooth);
    printf("BAD_VALUE   = %g\n", par->bad_value);
    printf("BKG_FRAC    = %g\n", par->bkg_frac);
    printf("THRESHOLD   = %g\n", par->threshold);
    printf("MIN_SCALE   = %g\n", par->min_scale);
    printf("MAX_SCALE   = %g\n", par->max_scale);
    printf("MIN_NGOOD   = %d\n", par->min_ngood);
    printf("N_REJ_LOW   = %d\n", par->nlow);
    printf("N_REJ_HIGH  = %d\n", par->nhigh);
    printf("MEDIAN      = %d\n", par->median);
    printf("VERBOSE     = %d\n", par->verbose);
    printf("-----------\n");
  }

// Check the parameters
  if (par->nsmooth > par->nbin)
  {
    printf("\n\tERROR: HIST_NBIN must not be smaller than NBIN_SMOOTH\n");
    exit(EXIT_FAILURE);
  }

  if (par->high < par->low)
  {
    printf("\n\tERROR: HIST_HIGH must not be smaller than HIST_LOW\n");
    exit(EXIT_FAILURE);
  }

  if (!(inpf=fopen(instrname, "r"))) errmess(instrname);

  while (!feof(inpf))
  {
    if (fgets(line, 200, inpf) == NULL) errmess(instrname);
    if (feof(inpf)) break;

    sscanf(line, "%s = %s", key, val);

    if (!strcasecmp(key, "END")) break;
    else if (!strcasecmp(key, "SAT_LEVEL")) sscanf(val, "%f", &par->sat_level);
    else if (!strcasecmp(key, "MIN_LEVEL")) sscanf(val, "%f", &par->min_level);
    else if (!strcasecmp(key, "NX"))        sscanf(val, "%d", &snx);
    else if (!strcasecmp(key, "NY"))        sscanf(val, "%d", &sny);
    else if (!strcasecmp(key, "MARG"))      sscanf(val, "%d", &smarg);
    else if (!strcasecmp(key, "FITS"))      continue;
    else if (!strcasecmp(key, "GAIN"))      continue;
    else printf("WARNING: Unknown parameter %s\n", key);
  }

  fclose(inpf);

  if (par->verbose > 1)
  {
    printf("Reading %s\n", instrname);
    printf("-----------\n");
    printf("SAT_LEVEL = %g\n", par->sat_level);
    printf("MIN_LEVEL = %g\n", par->min_level);
    printf("NX        = %d\n", snx);
    printf("NY        = %d\n", sny);
    printf("MARG      = %d\n", smarg);
    printf("-----------\n");
  }

  par->nx0=snx+smarg+smarg;
  par->ny0=sny+smarg+smarg;

  if (par->verbose > 1)
  {
    printf("NX0 =  %d\n", par->nx0);
    printf("NY0 =  %d\n", par->ny0);
    printf("-----------\n");
  }

  return;
}
/*** END ***/
