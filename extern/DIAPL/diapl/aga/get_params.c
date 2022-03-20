/*========================================================*/
/*                                                        */
/*  get_params.c        version 1.5.1   2005.07.08        */
/*  (aga)                                                 */
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

void get_params(char *parfname, char *instrname, PARAMS *par)
{
        char  line[201],
              key[200],
              val[200];
        int   snx, sny, smarg;
        FILE  *inpf;

  if (!(inpf = fopen(parfname, "r"))) errmess(parfname);

  while (!feof(inpf))
  {
    fgets(line, 200, inpf);
    if (feof(inpf)) break;

    sscanf(line, "%s = %s", key, val);

    if (!strcasecmp(key, "END")) break;
    else if (!strcasecmp(key, "KER_HW"))        sscanf(val, "%d", &par->kerhw);
    else if (!strcasecmp(key, "N_COMP"))        sscanf(val, "%d", &par->ncomp);
    else if (!strcasecmp(key, "DEG_1"))         sscanf(val, "%d", par->deg);
    else if (!strcasecmp(key, "SIG_GAUSS_1"))   sscanf(val, "%lf", par->sig);
    else if (!strcasecmp(key, "DEG_INC"))
      sscanf(val, "%d", &par->deg_inc);
    else if (!strcasecmp(key, "SIG_GAUSS_INC"))
      sscanf(val, "%lf", &par->sig_inc);
    else if (!strcasecmp(key, "BKG_DEG"))       sscanf(val, "%d", &par->bdeg);
    else if (!strcasecmp(key, "BAD_VALUE"))
      sscanf(val, "%lf", &par->bad_value);
    else if (!strcasecmp(key, "BAD_GROWRAD1"))
      sscanf(val, "%d", &par->bad_growrad1);
    else if (!strcasecmp(key, "BAD_GROWRAD2"))
      sscanf(val, "%d", &par->bad_growrad2);
    else if (!strcasecmp(key, "MIN_AREA"))
      sscanf(val, "%lf", &par->min_area);
    else if (!strcasecmp(key, "MAX_NITER"))     sscanf(val, "%d", &par->n_iter);
    else if (!strcasecmp(key, "N_SIG"))         sscanf(val, "%lf", &par->n_sig);
    else if (!strcasecmp(key, "MAX_CHI2"))
      sscanf(val, "%lf", &par->max_chi2);
    else if (!strcasecmp(key, "VERBOSE"))
      sscanf(val, "%d", &par->verbose);
    else if (!strcasecmp(key, "WDEG_SPATIAL"))  sscanf(val, "%d", &par->wdeg);
    else if (!strcasecmp(key, "DOM_HW"))        sscanf(val, "%d", &par->domhw);
    else if (!strcasecmp(key, "NDOM_X"))        sscanf(val, "%d", &par->ndom_x);
    else if (!strcasecmp(key, "NDOM_Y"))        sscanf(val, "%d", &par->ndom_y);
    else if (!strcasecmp(key, "DOM_THRESH"))
      sscanf(val, "%lf", &par->dom_thresh);
    else if (!strcasecmp(key, "N_SIG_DOM"))
      sscanf(val, "%lf", &par->n_sig_dom);
    else if (!strcasecmp(key, "DOMAIN_MODE"))
      sscanf(val, "%d", &par->domain_mode);
    else if (!strcasecmp(key, "MOHW"))          sscanf(val, "%d", &par->mohw);
    else if (!strcasecmp(key, "N_ITER_DOM"))
      sscanf(val, "%d", &par->n_iter_dom);
    else if (!strcasecmp(key, "MIN_AREA_DOM"))
      sscanf(val, "%lf", &par->min_area_dom);
    else if (!strcasecmp(key, "MIN_NKEEP"))
      sscanf(val, "%d", &par->min_nkeep);
    else printf("WARNING: Unknown parameter %s\n", key);
  }

  fclose(inpf);

  if (par->verbose > 1)
  {
    printf("Reading: %s\n", parfname);
    printf("-----------\n");
    printf("KER_HW        = %d\n", par->kerhw);
    printf("N_COMP        = %d\n", par->ncomp);
    printf("DEG_1         = %d\n", *par->deg);
    printf("SIG_GAUSS_1   = %g\n", *par->sig);
    printf("DEG_INC       = %d\n", par->deg_inc);
    printf("SIG_GAUSS_INC = %g\n", par->sig_inc);
    printf("BKG_DEG       = %d\n", par->bdeg);
    printf("BAD_VALUE     = %g\n", par->bad_value);
    printf("BAD_GROWRAD1  = %d\n", par->bad_growrad1);
    printf("BAD_GROWRAD2  = %d\n", par->bad_growrad2);
    printf("MIN_AREA      = %g\n", par->min_area);
    printf("MAX_NITER     = %d\n", par->n_iter);
    printf("N_SIG         = %g\n", par->n_sig);
    printf("MAX_CHI2      = %g\n", par->max_chi2);
    printf("WDEG_SPATIAL  = %d\n", par->wdeg);
    printf("DOM_HW        = %d\n", par->domhw);
    printf("NDOM_X        = %d\n", par->ndom_x);
    printf("NDOM_Y        = %d\n", par->ndom_y);
    printf("DOM_THRESH    = %g\n", par->dom_thresh);
    printf("N_SIG_DOM     = %g\n", par->n_sig_dom);
    printf("DOMAIN_MODE   = %d\n", par->domain_mode);
    printf("MOHW          = %d\n", par->mohw);
    printf("N_ITER_DOM    = %d\n", par->n_iter_dom);
    printf("MIN_AREA_DOM  = %g\n", par->min_area_dom);
    printf("MIN_NKEEP     = %d\n", par->min_nkeep);
    printf("VERBOSE       = %d\n", par->verbose);
    printf("-----------\n");
  }

  if (par->verbose > 1) printf("Reading %s\n", instrname);
  if (!(inpf=fopen(instrname, "r"))) errmess(instrname);

  while (!feof(inpf))
  {
    fgets(line, 200, inpf);
    if (feof(inpf)) break;

    sscanf(line, "%s = %s", key, val);

    if (!strcasecmp(key, "END")) break;
    else if (!strcasecmp(key, "SAT_LEVEL")) sscanf(val, "%lf", &par->sat_level);
    else if (!strcasecmp(key, "MIN_LEVEL")) sscanf(val, "%lf", &par->min_level);
    else if (!strcasecmp(key, "GAIN"))      sscanf(val, "%lf", &par->gain);
    else if (!strcasecmp(key, "NX"))        sscanf(val, "%d", &snx);
    else if (!strcasecmp(key, "NY"))        sscanf(val, "%d", &sny);
    else if (!strcasecmp(key, "MARG"))      sscanf(val, "%d", &smarg);
    else if (!strcasecmp(key, "FITS"))      sscanf(val, "%s", par->fits_ext);
    else printf("WARNING: Unknown parameter %s\n", key);
  }

  fclose(inpf);

  if (par->verbose > 1)
  {
    printf("-----------\n");
    printf("GAIN      = %g\n", par->gain);
    printf("SAT_LEVEL = %g\n", par->sat_level);
    printf("MIN_LEVEL = %g\n", par->min_level);
    printf("SNX       = %d\n", snx);
    printf("SNY       = %d\n", sny);
    printf("MARG      = %d\n", smarg);
    printf("FITS      = %s\n", par->fits_ext);
    printf("-----------\n");
  }

  par->nx0=snx+smarg+smarg;
  par->ny0=sny+smarg+smarg;

  par->nx=par->nx0-2*par->kerhw;
  par->ny=par->ny0-2*par->kerhw;

  if (par->fits_ext[0] == '"') strcpy(par->fits_ext, (par->fits_ext+1));
  if (par->fits_ext[strlen(par->fits_ext)-1] == '"')
    par->fits_ext[strlen(par->fits_ext)-1]='\0';

  if (par->verbose > 1)
  {
    printf("-----------\n");
    printf("NX=NX0-2*KER_HW= %d\n", par->nx);
    printf("NY=NY0-2*KER_HW= %d\n", par->ny);
    printf("NX0=SNX+2*MARG=  %d\n", par->nx0);
    printf("NY0=SNY+2*MARG=  %d\n", par->ny0);
    printf("-----------\n");
  }

  return;
}
/*** END ***/
