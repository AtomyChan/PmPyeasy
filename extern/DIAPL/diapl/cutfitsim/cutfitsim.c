/*========================================================*/
/*                                                        */
/*  cutfitsim.c         version 1.0.2   2005.05.06        */
/*                                                        */
/*  Copyright (C) 2005 by Wojtek Pych, CAMK PAN           */
/*  pych@camk.edu.pl                                      */
/*                                                        */
/*  Copyright (C) 2005 by Przemek Wozniak                 */
/*  wozniak@lanl.gov                                      */
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
#include <stdlib.h>
#include <string.h>

#include "errmess.h"
#include "swap.h" 
#include "pfitshead.h"

#define verbose 0

/*--------------------------------------------------------*/
char *strip_path(char *fname)
{
        char *cp;

  cp=strrchr(fname, '/');

  if (cp) return (cp+1);
  else    return fname;
}
/*--------------------------------------------------------*/
short from_header(int hsize, char **header, int *naxis1, int *naxis2,
                  double *Fbzero, double *Fbscale)
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
void change_n(int hsize, char **header, char *keyword, int naxis)
{
        char  card[CARD_SIZE]; 
        int   p;
        
  p=get_FITS_key(hsize, header, keyword, card);

  sprintf(card, "%-8s= %20d /", keyword, naxis);
  memcpy(header[p], card, strlen(card));

  return;
}
/*--------------------------------------------------------*/
int modify_header(int hsize, char ***header, int n11, int n12, int n21, int n22,
                  char *iname)
{
        char  value[CARD_SIZE],
              tmp[RECORD_SIZE],
              **loc_head;
        int   p,
              i;

  loc_head=*header;

  change_n(hsize, loc_head, "NAXIS1", n12-n11+1);
  change_n(hsize, loc_head, "NAXIS2", n22-n21+1);

  if ((p=get_FITS_key(hsize, loc_head, "END", value)) == -1)
  {
    printf("ERROR! modify_header(END)\n");
    return(-1);
  } 

  memset(tmp, ' ', CARD_SIZE);
  sprintf(tmp, "HISTORY   subframe [%d:%d,%d:%d] of %s",
          n11, n12, n21, n22, iname);
  tmp[strlen(tmp)]=' ';
  memcpy(value, tmp, CARD_SIZE);

  if (!((p+1)%RECORD_CARDS))
  {
    hsize+=RECORD_CARDS;
    if (!(loc_head=(char **)realloc(loc_head, hsize*sizeof(char *))))
      errmess("realloc(loc_head)");
    for (i=1; i<=RECORD_CARDS; i++)
    {
      if (!(loc_head[hsize-i]=(char *)malloc(CARD_SIZE)))
        errmess("malloc(loc_head[hsize-1])");
      memset(loc_head[hsize-i], ' ', CARD_SIZE);
    }
  }

  memcpy(loc_head[p+1], loc_head[p], CARD_SIZE);
  memcpy(loc_head[p], value, CARD_SIZE);

  *header=loc_head;
 
  return(hsize);
}
/*--------------------------------------------------------*/
int main(int argc, char *argv[])
{
        void    *pixbuf;
        char    *iname,          // input FITS file name
                *oname,          // output FITS file name
                **header,        // FITS header
                cimzero,
                buf[2*RECORD_SIZE],
                *cpixbuf;
        short   bitpix,
                simzero,
                *spixbuf;
        int     bytepix,
                ns,             // swap needed when != 0
                naxis1, naxis2,
                pixbufsize,
                hsize,
                n11, n12, n21, n22, nx, ny, i, j, k,
                dn11, dn12, dn21, dn22;
        long    ofs,
                limzero,
                *lpixbuf;
        double   Fbzero,
                Fbscale,
                imzero,
                *fpixbuf;
        double  *dpixbuf;
        FILE    *inf, *outf;

  if (argc != 7)
  {
    printf("\n\tUSAGE: %s infile outfile n11 n12 n21 n22\n", argv[0]);
    exit(1);
  }

  iname=argv[1];
  oname=argv[2];

  if (sscanf(argv[3], "%d", &n11) != 1)
  { printf("ERROR! in 'n11' argument\n"); exit(2);  }

  if (sscanf(argv[4], "%d", &n12) != 1)
  { printf("ERROR! in 'n12' argument\n"); exit(3);  }

  if (sscanf(argv[5], "%d", &n21) != 1)
  { printf("ERROR! in 'n21' argument\n"); exit(4);  }

  if (sscanf(argv[6], "%d", &n22) != 1)
  { printf("ERROR! in 'n22' argument\n"); exit(5);  }

  if ((n11 >= n12)  ||  (n21 >= n22))
  {
    printf("ERROR! bad indices:  n11 >= n12  or  n21 >= n22 !\n");
    exit(6);
  }

  if (verbose)
  {
    printf("Input file: %s\n", iname);
    printf("Output file: %s\n", oname);
    printf("x1= %d   x2= %d   y1= %d   y2= %d\n", n11, n12, n21, n22);
  }

  if (!(inf=fopen(iname, "r")))   errmess(iname);

  hsize=read_FITS_header(inf, &header);
  ofs=(long)hsize*CARD_SIZE;

  bitpix=from_header(hsize, header, &naxis1, &naxis2, &Fbzero, &Fbscale);
  if (bitpix == -1)
  {
    printf("ERROR! Reading header\n");
    exit(7);
  }
  bytepix=abs(bitpix)/8;

  if (verbose)
  {
    printf("%s: %d x %d\n", iname, naxis1, naxis2);
    printf("header size= %d   offset= %ld\n", hsize, ofs);
  }

  hsize=modify_header(hsize, &header, n11, n12, n21, n22, iname);

  if (!(outf=fopen(oname, "w")))  errmess(oname);

  write_FITS_header(outf, hsize, header);

/* fancy reading for all sorts of cuts colliding with boundaries */

  if (n21 < 1)  dn21 = -n21+1;
  else          dn21 = 0;

  if (n22 > naxis2) dn22 = n22 - naxis2;
  else              dn22 = 0;

  if (n11 < 1)  dn11 = -n11+1;
  else          dn11 = 0;

  if (n12 > naxis1) dn12 = n12 - naxis1;
  else              dn12 = 0;

  nx=n12-n11+1;
  ny=n22-n21+1;

  if ((n11 > naxis1) || (n12 < 1) || (n21 > naxis2) || (n22 < 1))
  {
    dn21 = ny;
    dn22 = 0;
  }

  if (verbose)
  {
    printf("dn11= %d   dn12= %d   dn21= %d   dn22= %d\n",
            dn11, dn12, dn21, dn22);
    printf("nx= %d   ny= %d\n", nx, ny);
  }

  pixbufsize=(long)nx*bytepix;
  imzero=-Fbzero/Fbscale;

  ns=needswap();

  cpixbuf = NULL;
  spixbuf = NULL;
  lpixbuf = NULL;
  fpixbuf = NULL;
  dpixbuf = NULL;

  switch(bitpix)
  {
    case   8:
      if (!(cpixbuf=(char *)malloc(pixbufsize)))
        errmess("calloc(cpixbuf)");
      pixbuf=(void *)cpixbuf;
      cimzero=(char)imzero;
      if (ns) swap_bytes(&cimzero, 1, bytepix);
      for (i=0; i<nx; i++) cpixbuf[i]=cimzero;
      break;
    case  16:
      if (!(spixbuf=(short *)malloc(pixbufsize)))
        errmess("calloc(spixbuf)"); 
      pixbuf=(void *)spixbuf;
      simzero=(short)imzero;
      if (ns) swap_bytes(&simzero, 1, bytepix);
      for (i=0; i<nx; i++) spixbuf[i]=simzero;
      break;
    case  32:
      if (!(lpixbuf=(long *)malloc(pixbufsize)))
        errmess("calloc(lpixbuf)"); 
      pixbuf=(void *)lpixbuf;
      limzero=(long)imzero;
      if (ns) swap_bytes(&limzero, 1, bytepix);
      for (i=0; i<nx; i++) lpixbuf[i]=limzero;
      break;
    case -32:
      if (!(fpixbuf=(double *)malloc(pixbufsize)))
        errmess("calloc(fpixbuf)"); 
      pixbuf=(void *)fpixbuf;
      if (ns) swap_bytes(&imzero, 1, bytepix);
      for (i=0; i<nx; i++) fpixbuf[i]=imzero;
      break;
    case -64:
      if (!(dpixbuf=(double *)malloc(pixbufsize)))
        errmess("calloc(dpixbuf)"); 
      pixbuf=(void *)dpixbuf;
      if (ns) swap_bytes(&imzero, 1, bytepix);
      for (i=0; i<nx; i++) dpixbuf[i]=imzero;
      break;
    default:
      printf("BITPIX not conforming FITS standard - not supported\n");
      exit(8);
  }

  for (j=1; j<=dn21; j++ )
  {
    if (verbose) printf("Zeros: j= %d\n", j);

    if (fwrite(pixbuf, bytepix, nx, outf) != nx)
    {
      printf("ERROR! writing pixel data to %s\nExiting...\n", oname);
      exit(9);
    }
  }

  ofs+=(long)(n11+dn11-1)*bytepix;
  if (verbose) printf("n11+dn11-1= %d\n", n11+dn11-1);
  for (j=n21+dn21; j<=n22-dn22; j++)
  {
    if (verbose) printf("j= %d   ofs= %ld\n", j, ofs);

    fseek(inf, ofs+(long)(j-1)*naxis1*bytepix, SEEK_SET);

    if (fread(pixbuf+bytepix*dn11, bytepix, nx-dn11-dn12, inf) != nx-dn11-dn12)
    {
      printf("ERROR! reading pixel data from %s\nExiting...\n", iname);
      exit(10);
    }

    if (fwrite(pixbuf, bytepix, nx, outf) != nx)
    {
      printf("ERROR! writing pixel data to %s\nExiting...\n", oname);
      exit(11);
    }
  }

  fclose(inf);

  switch(bitpix)
  {
    case   8: for (i=0; i<nx; i++) cpixbuf[i]=cimzero;
      break;
    case  16: for (i=0; i<nx; i++) spixbuf[i]=simzero;
      break;
    case  32: for (i=0; i<nx; i++) lpixbuf[i]=limzero;
      break;
    case -32: for (i=0; i<nx; i++) fpixbuf[i]=imzero;
      break;
    case -64: for (i=0; i<nx; i++) dpixbuf[i]=imzero;
      break;
    default:  printf("BITPIX - this should not happen\n");
      exit(12);
  }

  for (j=1; j<=dn22; j++)
  {
    if (fwrite(pixbuf, bytepix, nx, outf) != nx)
    {
      printf("ERROR! writing pixel data to %s\nExiting...\n", oname);
      exit(13);
    }
  }

  free(pixbuf);

  k=((long)(nx)*(long)(ny)*bytepix)%RECORD_SIZE;
  k=RECORD_SIZE-k;
  if (k != RECORD_SIZE)
  {
    memset(buf, 0, k);
    if (fwrite(buf, 1, k, outf) != k)
    {
      printf("ERROR! writing pixel data to %s\nExiting...\n", oname);
      exit(14);
    }
  }

  fclose(outf);

  return(0);
}
/*** END ***/
