/*========================================================*/
/*                                                        */
/*  get_params.c        version 1.4.0   2005.05.06        */
/*  (getvar)                                              */
/*                                                        */
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
#include <string.h>

#include "defs.h"
#include "errmess.h"

void get_params(char *parfname, char *instrname, PAR_STRUCT *par)
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
    else if (!strcasecmp(key, "APRAD"))     sscanf(val, "%f", &par->aprad);
    else if (!strcasecmp(key, "FITRAD"))    sscanf(val, "%f", &par->fitrad);
    else if (!strcasecmp(key, "NORMRAD"))   sscanf(val, "%f", &par->normrad);
    else if (!strcasecmp(key, "ANRAD1"))    sscanf(val, "%f", &par->anrad1);
    else if (!strcasecmp(key, "ANRAD2"))    sscanf(val, "%f", &par->anrad2);
    else if (!strcasecmp(key, "BKG_MODE"))  sscanf(val, "%d", &par->bkg_mode);
    else if (!strcasecmp(key, "NSIG_BKG"))  sscanf(val, "%f", &par->nsig_bkg);
    else if (!strcasecmp(key, "BAD_VALUE")) sscanf(val, "%f", &par->bad_value);
    else if (!strcasecmp(key, "ERR_CODE"))  sscanf(val, "%f", &par->err);
    else if (!strcasecmp(key, "SMHW"))      sscanf(val, "%d", &par->smhw);
    else if (!strcasecmp(key, "MOHW"))      sscanf(val, "%d", &par->mohw);
    else if (!strcasecmp(key, "C_MIN"))     sscanf(val, "%f", &par->c_min);
    else if (!strcasecmp(key, "NOBJ_INIT")) sscanf(val, "%d", &par->nobj_init);
    else if (!strcasecmp(key, "ID_RAD"))    sscanf(val, "%f", &par->id_rad);
    else if (!strcasecmp(key, "NCONS_VAR1")) sscanf(val, "%d", &par->ncons_var1);
    else if (!strcasecmp(key, "NPTS_VAR2")) sscanf(val, "%d", &par->npts_var2);
    else if (!strcasecmp(key, "NSIG_VAR1")) sscanf(val, "%f", &par->nsig_var1);
    else if (!strcasecmp(key, "NSIG_VAR2")) sscanf(val, "%f", &par->nsig_var2);
    else if (!strcasecmp(key, "LIM_RATIO")) sscanf(val, "%f", &par->lim_ratio);
    else if (!strcasecmp(key, "ISOHW"))     sscanf(val, "%d", &par->isohw);
    else if (!strcasecmp(key, "ISO_OFF"))   sscanf(val, "%f", &par->iso_off);
    else if (!strcasecmp(key, "ISO_SLO"))   sscanf(val, "%f", &par->iso_slo);
    else if (!strcasecmp(key, "BAD_MARGIN")) sscanf(val, "%f", &par->bad_margin);
    else if (!strcasecmp(key, "CENTER_HW")) sscanf(val, "%d", &par->center_hw);
    else if (!strcasecmp(key, "FILTER_HW")) sscanf(val, "%d", &par->filter_hw);
    else if (!strcasecmp(key, "FWHM_FRAC")) sscanf(val, "%f", &par->fwhm_frac);
    else if (!strcasecmp(key, "VERBOSE"))   sscanf(val, "%d", &par->verbose);
    else printf("WARNING: Unknown parameter %s\n", key);
  }

  fclose(inpf);

  if (par->verbose > 2)
  {
    printf("Reading %s:\n", parfname);
    printf("-----------\n");
    printf("APRAD       = %g\n", par->aprad);
    printf("FITRAD      = %g\n", par->fitrad);
    printf("NORMRAD     = %g\n", par->normrad);
    printf("ANRAD1      = %g\n", par->anrad1);
    printf("ANRAD2      = %g\n", par->anrad2);
    printf("BKG_MODE    = %d\n", par->bkg_mode);
    printf("NSIG_BKG    = %g\n", par->nsig_bkg);
    printf("BAD_VALUE   = %g\n", par->bad_value);
    printf("ERR_CODE    = %g\n", par->err);
    printf("SMHW        = %d\n", par->smhw);
    printf("MOHW        = %d\n", par->mohw);
    printf("C_MIN       = %g\n", par->c_min);
    printf("NOBJ_INIT   = %d\n", par->nobj_init);
    printf("ID_RAD      = %g\n", par->id_rad);
    printf("NCONS_VAR1  = %d\n", par->ncons_var1);
    printf("NPTS_VAR2   = %d\n", par->npts_var2);
    printf("NSIG_VAR1   = %g\n", par->nsig_var1);
    printf("NSIG_VAR2   = %g\n", par->nsig_var2);
    printf("LIM_RATIO   = %g\n", par->lim_ratio);
    printf("ISOHW       = %d\n", par->isohw);
    printf("ISO_OFF     = %g\n", par->iso_off);
    printf("ISO_SLO     = %g\n", par->iso_slo);
    printf("BAD_MARGIN  = %g\n", par->bad_margin);
    printf("CENTER_HW   = %d\n", par->center_hw);
    printf("FILTER_HW   = %d\n", par->filter_hw);
    printf("FWHM_FRAC   = %g\n", par->fwhm_frac);
    printf("VERBOSE     = %d\n", par->verbose);
    printf("-----------\n");
  }

  if (!(inpf=fopen(instrname, "r"))) errmess(instrname);

  while (!feof(inpf))
  {
    fgets(line, 200, inpf);
    if (feof(inpf)) break;

    sscanf(line, "%s = %s", key, val);

    if (!strcasecmp(key, "END")) break;
    else if (!strcasecmp(key, "SAT_LEVEL")) sscanf(val, "%f", &par->sat_level);
    else if (!strcasecmp(key, "MIN_LEVEL")) sscanf(val, "%f", &par->min_level);
    else if (!strcasecmp(key, "GAIN"))      sscanf(val, "%f", &par->gain);
    else if (!strcasecmp(key, "NX"))        continue;
    else if (!strcasecmp(key, "NY"))        continue;
    else if (!strcasecmp(key, "MARG"))      continue;
    else printf("WARNING: Unknown parameter %s\n", key);
  }

  fclose(inpf);

  if (par->verbose)
  {
    printf("Reading %s\n", instrname);
    printf("-----------\n");
    printf("GAIN      = %g\n", par->gain);
    printf("SAT_LEVEL = %g\n", par->sat_level);
    printf("MIN_LEVEL = %g\n", par->min_level);
    printf("-----------\n");
  }

  return;
}
/*** END ***/
