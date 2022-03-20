/*========================================================*/
/*                                                        */
/*  get_params.c        version 1.9.0   2005.05.06        */
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

#include <ctype.h>
#include <stdio.h>
#include <string.h>

#include "defs.h"
#include "errmess.h"

void get_params(char *parfname, char *instrname, PAR_STRUCT *par)
{
        char  line[201],
              key[200],
              val[200];
        FILE  *inpf;

  if (!(inpf = fopen(parfname, "r"))) errmess(parfname);

  while (!feof(inpf))
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
    else if (!strcasecmp(key, "DBFORM"))    sscanf(val, "%c", &par->dbf);
    else if (!strcasecmp(key, "VERBOSE"))   sscanf(val, "%d", &par->verbose);
    else printf("WARNING: Unknown parameter %s\n", key);
  }

  fclose(inpf);

  par->dbf=toupper(par->dbf);

  if (par->verbose > 1)
  {
    printf("APRAD    = %g\n", par->aprad);
    printf("FITRAD   = %g\n", par->fitrad);
    printf("NORMRAD  = %g\n", par->normrad);
    printf("ANRAD1   = %g\n", par->anrad1);
    printf("ANRAD2   = %g\n", par->anrad2);
    printf("BKG_MODE = %d\n", par->bkg_mode);
    printf("NSIG_BKG = %g\n", par->nsig_bkg);
    printf("BAD_VALUE= %g\n", par->bad_value);
    printf("ERR_CODE = %g\n", par->err);
    printf("DBFORM   = %c\n", par->dbf);
    printf("VERBOSE  = %d\n", par->verbose);
    printf("---------\n");
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
