/*========================================================*/
/*                                                        */
/*  mkushortone.c       version 1.0.3     2006.10.24      */
/*                                                        */
/*  Copyright (C) 2006 by Wojtek Pych, CAMK PAN           */
/*  pych@camk.edu.pl                                      */
/*                                                        */
/*  Written for GNU project C and C++ Compiler            */
/*                                                        */
/*  Create 16-bit square FITS image with values 1         */
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

#include "errmess.h"
#include "pfitsio.h"

void usage(void);

/*--------------------------------------------------------*/
void usage()
{
  printf("\n\t USAGE: mkushortone [n output_file]\n\n");
  printf("  n - size of the image in pixels\n\n");
  exit(1);
}
/*--------------------------------------------------------*/
int create_header(int n, char ***header)
{
        char  hline[CARD_SIZE],
              **lheader;
        int   i;

  if (!(lheader=(char **)calloc(6, sizeof(char *)))) errmess("calloc(lheader)");
  for (i=0; i<6; i++)
  {
    if (!(lheader[i]=(char *)calloc(CARD_SIZE, sizeof(char))))
      errmess("calloc(lheader[i])");
    memset(lheader[i], ' ', CARD_SIZE);
  }

  sprintf(hline, "SIMPLE  =                    T");
  memcpy(lheader[0], hline, strlen(hline));

  sprintf(hline, "BITPIX  =                   16 /");
  memcpy(lheader[1], hline, strlen(hline));

  sprintf(hline, "NAXIS   =                    2 /");
  memcpy(lheader[2], hline, strlen(hline));

  sprintf(hline, "NAXIS1  =                 %4d /", n);
  memcpy(lheader[3], hline, strlen(hline));

  sprintf(hline, "NAXIS2  =                 %4d /", n);
  memcpy(lheader[4], hline, strlen(hline));

  sprintf(hline, "END");
  memcpy(lheader[5], hline, strlen(hline));

  *header=lheader;

  return(6);
}
/*--------------------------------------------------------*/
int main(int argc, char *argv[])
{
  const int   bytepix=2;
        char  *oname,
              **header,
              *tmp;
        short *row;
        int   i,
              nx,
              hsize,
              ls;
        long  size;
        FILE  *outf;

  if (argc == 3)
  {
    nx=atoi(argv[1]);
    oname=argv[2];
  }
  else usage();

  hsize=create_header(nx, &header);

  if (!(row=(short *)calloc((size_t)nx, sizeof(short)))) errmess("calloc(row)");
  for (i=0; i<nx; i++) row[i]=1;
  if (needswap()) swap_bytes(row, (size_t)nx, bytepix);

  if (!(outf=fopen(oname, "w")))  errmess(oname);

  write_FITS_header(outf, hsize, header);
  for (i=0; i<nx; i++)
  {
    if (fwrite(row, (size_t)bytepix, (size_t)nx, outf) != (size_t)nx)
    {
      printf("ERROR! writing image to %s failed\n", oname);
      exit(EXIT_FAILURE);
    }
  }

  size=nx*nx*bytepix;
  if ((ls=size%RECORD_SIZE) != 0) ls=RECORD_SIZE-ls;
  if (!(tmp=(char *)malloc((size_t)ls))) errmess("malloc(tmp)");
  memset(tmp, 0, (size_t)ls);
  if (fwrite(tmp, 1, (size_t)ls, outf) != (size_t)ls)
  {
    printf("ERROR! writing image to %s failed\n", oname);
    exit(EXIT_FAILURE);
  }
  fclose(outf);

  for (i=0; i<hsize; i++) free(header[i]);
  free(header);
  free(row);
  free(tmp);

  return(0);
}
/*** END ***/
