/*========================================================*/
/*                                                        */
/*  init_difimages.c    version 1.4.0   2005.04.19        */
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
#include <string.h>

#include "defs.h"
#include "funcs.h"
#include "errmess.h"
#include "pfitsio1.h"

void init_difimages(char **inpfname, char **outfname, long *headlen,
                    int nim, PARAMS par)
{
        char  **header,
              value[VALUE_SIZE];
        int   i, j, k, bpx, nx, ny,
              hsize;
        double *buf, tmp;
        FILE  *fp;

/* make a buffer full of bad values */
  if (!(buf = (double *)malloc(par.nx0*sizeof(double)))) errmess("malloc(buf)");

  tmp = par.bad_value;
  if (needswap()) swap_bytes(&tmp, 1, sizeof(double));
  for (i=0; i<par.nx0; i++) buf[i] = tmp;

  for (i=0; i<nim; i++)
  {
    if (!(fp=fopen(inpfname[i], "r"))) errmess(inpfname[i]);
    hsize=read_FITS_header(fp, &header);
    fclose(fp);

    if (hsize * CARD_SIZE % RECORD_SIZE)
      headlen[i]=((hsize * CARD_SIZE / RECORD_SIZE)+1)*RECORD_SIZE;
    else
      headlen[i]=hsize * CARD_SIZE;

    if (get_FITS_key(hsize, header, "NAXIS1", value) == -1)
      errmess("Cannot find NAXIS1 in a FITS header");
    sscanf(value, "%d", &nx);

    if (get_FITS_key(hsize, header, "NAXIS2", value) == -1)
      errmess("Cannot find NAXIS2 in a FITS header");
    sscanf(value, "%d", &ny);

    if (get_FITS_key(hsize, header, "BITPIX", value) == -1)
      errmess("Cannot find BITPIX in a FITS header");
    sscanf(value, "%d", &bpx);

    if (nx != par.nx0)
    {
      printf("ERROR: aga-init_difimages(%s): ", inpfname[i]);
      printf("NAXIS1: expected= %d    got= %d\n", par.nx0, nx);
      exit(2);
    }

    if (ny != par.ny0)
    {
      printf("ERROR: aga-init_difimages(%s): ", inpfname[i]);
      printf("NAXIS2: expected= %d    got %d\n", par.ny0, ny);
      exit(3);
    }

/****
    if (bpx != -32)
    {
      printf("init_difimages: %s error: ", inpfname[i]);
      printf("BITPIX must be -32\n");
      exit(4);
    }
****/

/* modify header */
    if ((j=get_FITS_key(hsize, header, "BITPIX", value)) == -1)
    {
      printf("init_difimages: BITPIX not found in the FITS header\n");
      exit(2);
    }
    strcpy(header[j], "BITPIX  =                  -32");
    strcat(header[j], " / Bits per pixel                                 ");

    if (!(fp=fopen(outfname[i], "w"))) errmess(outfname[i]);
    write_FITS_header(fp, hsize, header);

    for (k=0; k<par.ny0; k++)
      fwrite(buf, sizeof(double), par.nx0, fp);

    fclose(fp);

    for (k=0; k<hsize; k++) free(header[k]);
    free(header);
  }

  free(buf);

  return;
}
/*** END ***/
