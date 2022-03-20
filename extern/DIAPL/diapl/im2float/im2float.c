/*========================================================*/
/*                                                        */
/*  im2float.c          version 0.1.1   2005.07.08        */
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
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "errmess.h"
#include "pfitsio.h"

/*--------------------------------------------------------*/
short from_header(int hsize, char **header, int *naxis1, int *naxis2,
                  float *Fbzero, float *Fbscale)
{
        char  value[VALUE_SIZE],
              val[VALUE_SIZE];
        short naxis,
              lbitpix;
        int   p;

  if ((p=get_FITS_key(hsize, header, "SIMPLE", value)) == -1)
  {
    printf("ERROR! SIMPLE not found in the header\n");
    return(-1);
  }
  if (p != 0)
    printf("WARNING! Header does not conform the FITS standard (SIMPLE)\n");
  sscanf(value, "%s", val);
  if ((strcmp(val, "T")) && (strncmp(val, "T/", 2)))
    printf("WARNING! File does not conform the FITS standard (SIMPLE)\n");

  if ((p=get_FITS_key(hsize, header, "BITPIX", val)) == -1)
  {
    printf("ERROR! BITPIX not found in the header\n");
    return(-1);
  }
  if (p != 1)
    printf("WARNING! header does not conform to FITS standard (BITPIX)\n");
  sscanf(val, "%hd", &lbitpix);

  if ((p=get_FITS_key(hsize, header, "NAXIS", val)) == -1)
  {
    printf("ERROR! NAXIS not found in the header\n");
    return(-1);
  }
  if (p != 2)
    printf("WARNING! header does not conform to FITS standard (NAXIS)\n");

  sscanf(val, "%hd", &naxis);
  if (naxis != 2)
  {
    printf("ERROR! Image must be 2-dimensional\n");
    return(-1);
  }

  if ((p=get_FITS_key(hsize, header, "NAXIS1", val)) == -1)
  {
    printf("ERROR! NAXIS1 not found in the header\n");
    return(-1);
  }
  if (p != 3)
    printf("WARNING! header does not conform to FITS standard (NAXIS1)\n");

  sscanf(val, "%d", naxis1);

  if ((p=get_FITS_key(hsize, header, "NAXIS2", val)) == -1)
  {
    printf("ERROR! NAXIS2 not found in the header\n");
    return(-1);
  }
  if (p != 4)
    printf("WARNING! header does not conform to FITS standard (NAXIS2)\n");

  sscanf(val, "%d", naxis2);

  if ((p=get_FITS_key(hsize, header, "BZERO", val)) != -1)
    sscanf(val, "%lf", Fbzero);
  else
    *Fbzero=0.0;
  if ((p=get_FITS_key(hsize, header, "BSCALE", val)) != -1)
    sscanf(val, "%lf", Fbscale);
  else
    *Fbscale=1.0;

  return(lbitpix);
}
/*--------------------------------------------------------*/
int modify_header(int hsize, char ***header)
{
        char  value[CARD_SIZE],
              **loc_head;
        int   p;

  loc_head=*header;

  if ((p=get_FITS_key(hsize, loc_head, "BITPIX", value)) == -1)
  {
    printf("ERROR! modify_header(): BITPIX not found in the header\n");
    return(-1);
  }

  memcpy(loc_head[p], "BITPIX  =                  -32", 30);

  hsize=del_header_card(hsize, &loc_head, "BZERO");
  hsize=del_header_card(hsize, &loc_head, "BSCALE");

  *header=loc_head;
 
  return(hsize);
}
/*--------------------------------------------------------*/
int main(int argc, char *argv[])
{
        char  *iname,
              *oname,
              **header;
        short bitpix;
        int   hsize,
              naxis1, naxis2,
              i;
        float Fbzero,
              Fbscale,
              **im;
        FILE  *iof;

  if (argc != 3)
  {
    printf("\n\tUSAGE: im2float  infile  outfile\n\n");
    exit(1);
  }

  iname=argv[1];
  oname=argv[2];

  if (!(iof=fopen(iname, "r"))) errmess(iname);
  hsize=read_FITS_header(iof, &header);
  fclose(iof);

  bitpix=from_header(hsize, header, &naxis1, &naxis2, &Fbzero, &Fbscale);
  if (bitpix == -1)
  {
    printf("ERROR! Reading header\n");
    exit(7);
  }

  for (i=0; i<hsize; i++) free(header[i]);
  free(header);

  if (bitpix == -32)
  {
    if (symlink(iname, oname)) errmess(oname);
    return(0);
  }

  im=read_FITS_2Dfile(iname, 's', &hsize, &header, &naxis1, &naxis2);

  if ((hsize=modify_header(hsize, &header)) == -1) return(-1);

  write_FITS_2Dfile(oname, hsize, header, naxis1, naxis2, 4, (void **)im);

  for (i=0; i<naxis2; i++) free(im[i]);
  free(im);

  for (i=0; i<hsize; i++) free(header[i]);
  free(header);

  return(0);
}
/*** END ***/
