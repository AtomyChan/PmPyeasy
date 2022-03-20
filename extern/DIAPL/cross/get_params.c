/*========================================================*/
/*                                                        */
/*  get_params.c        version 1.3.6   2005.05.12        */
/*  (cross)                                               */
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

void get_params(char *fname, char *instrname, PARAMS *par)
{
        char  line[201],
              key[200],
              val[200];
        int   i;
        FILE *inf;

  if (!(inf=fopen(fname, "r"))) errmess(fname);

  for (i=0; !feof(inf); i++)
  {
    fgets(line, 200, inf);
    if (feof(inf)) break;

    sscanf(line, "%s = %s", key, val);

    if (!strcasecmp(key, "END")) break;
    else if (!strcasecmp(key, "X_NBIN"))    sscanf(val, "%d", &par->x_nbin);
    else if (!strcasecmp(key, "Y_NBIN"))    sscanf(val, "%d", &par->y_nbin);
    else if (!strcasecmp(key, "NX_CUT"))    sscanf(val, "%d", &par->nx_cut);
    else if (!strcasecmp(key, "NY_CUT"))    sscanf(val, "%d", &par->ny_cut);
    else if (!strcasecmp(key, "BKG_FRAC"))  sscanf(val, "%f", &par->bkg_frac);
    else if (!strcasecmp(key, "VERBOSE"))   sscanf(val, "%d", &par->verbose);
    else printf("WARNING: Unknown parameter %s\n", key);
  }

  fclose(inf);

  if ((par->nx_cut%2) || (par->ny_cut%2))
  {
    printf("ERROR! get_params(): NX_CUT and NY_CUT must be a power of 2\n");
    exit(3);
  }

  if (par->verbose)
  {
    printf("Reading %s\n", fname);
    printf("-----------\n");
    printf("X_NBIN  = %d\n", par->x_nbin);
    printf("Y_NBIN  = %d\n", par->y_nbin);
    printf("NX_CUT  = %d\n", par->nx_cut);
    printf("NY_CUT  = %d\n", par->ny_cut);
    printf("BKG_FRAC= %g\n", par->bkg_frac);
    printf("VERBOSE = %d\n", par->verbose);
    printf("-----------\n");
  }

  if (!(inf=fopen(instrname, "r"))) errmess(instrname);

  for (i=0; !feof(inf); i++)
  {
    fgets(line, 200, inf);
    if (feof(inf)) break;

    sscanf(line, "%s = %s", key, val);

    if (!strcasecmp(key, "END")) break;
    else if (!strcasecmp(key, "SAT_LEVEL")) sscanf(val, "%f", &par->sat_level);
    else if (!strcasecmp(key, "MIN_LEVEL")) sscanf(val, "%f", &par->min_level);
    else if (!strcasecmp(key, "GAIN"))      continue;
    else if (!strcasecmp(key, "NX"))        continue;
    else if (!strcasecmp(key, "NY"))        continue;
    else if (!strcasecmp(key, "MARG"))      continue;
    else if (!strcasecmp(key, "FITS"))      continue;
    else printf("WARNING: Unknown parameter %s\n", key);
  }

  fclose(inf);

  if (par->verbose)
  {
    printf("Reading %s\n", instrname);
    printf("-----------\n");
    printf("SAT_LEVEL = %g\n", par->sat_level);
    printf("MIN_LEVEL = %g\n", par->min_level);
    printf("-----------\n");
  }

  return;
}
/*** END ***/
