/*========================================================*/
/*                                                        */
/*  pfitsio1.c      version 2.6.0         2006.10.25      */
/*                                                        */
/*  Copyright (C) 2006 by Wojtek Pych, CAMK PAN           */
/*                                                        */
/*  Written for GNU project C and C++ Compiler            */
/*  Private FITS I/O library.                             */
/*                                                        */
/*========================================================*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "pfitsio1.h"
#include "errmess.h"

double nearbyint(double x); // there may be a bug in current version of math.h

/*--------------------------------------------------------*/
/* Returns 1 on Intel architecture machines, 0 on Suns    */
/* Original function from DIA package by P. Wozniak       */
/*--------------------------------------------------------*/
int needswap(void)
{
        union { short s; char c[2]; } u;

  u.s=1;
  return((int)(u.c[0]));
}
/*--------------------------------------------------------*/
void swap_bytes(void *ptr, size_t ndata, int nbytes)
{
        char    c,
                *data;
        int     j, k;
        size_t  i, m;

  data=(char *)ptr;
  ndata*=nbytes;
  for (i=0; i<ndata; i+=nbytes)
  {
    k=nbytes/2;
    for (j=0; j<k; j++)
    {
      m=i+nbytes-1-j;
      c=data[i+j];
      data[i+j]=data[m];
      data[m]=c;
    }
  }

  return;
}
/*--------------------------------------------------------*/
/*  Read header card value with comment                   */
/*  Return card number (starting from 0)                  */
/*  or -1 on error.                                       */
/*--------------------------------------------------------*/
int get_FITS_key(int ncards, char **header, char *keyword, char *value)
{
        char  lkeyword[KEYWORD_SIZE];
        int   i;              /*  loop numerator          */

  if (strlen(keyword) > KEYWORD_SIZE)
  {
    printf("WARNING! get_FITS_key(%s): invalid keyword\n", keyword);
    return(-1);
  }

  for (i=0; i<(int)strlen(keyword); i++)            lkeyword[i]=keyword[i];
  for (i=(int)strlen(keyword); i<KEYWORD_SIZE; i++) lkeyword[i]=' ';

  for (i=0; i<ncards; i++)
  {
    if (!strncmp(header[i], lkeyword, KEYWORD_SIZE))
    {
      memcpy(value, header[i]+10, VALUE_SIZE);  // value may start in column 11
      break;
    }
  }

  if (i == ncards) return(-1);

  return(i);
}
/*--------------------------------------------------------*/
/*  Return: image size in bytes (without the filling)     */
/*          0 no image                                    */
/*         -1 error                                       */
/*--------------------------------------------------------*/
size_t FITS_image_size(int ncards, char **header)
{
        char    val[VALUE_SIZE],
                keyword[KEYWORD_SIZE];
        int     i,
                p,
                bitpix,
                bytepix,
                naxis,
                naxisn;
        size_t  size;

  if ((p=get_FITS_key(ncards, header, "BITPIX", val)) == -1)
  {
    printf("\n\tERROR! FITS_image_size(): BITPIX not found in the header\n");
    return(-1);
  }
  if (p != 1)
    printf("WARNING! header does not conform to FITS standard (BITPIX)\n");
  sscanf(val, "%d", &bitpix);
  bytepix=abs(bitpix)/8;

  strcpy(keyword, "NAXIS");
  if ((p=get_FITS_key(ncards, header, keyword, val)) == -1)
  {
    printf("\n\tERROR! FITS_image_size(): NAXIS not found in the header\n");
    return(-1);
  }
  if (p != 2)
    printf("WARNING! header does not conform to FITS standard (NAXIS)\n");

  sscanf(val, "%d", &naxis);
  if (naxis == 0) return(0);

  size=(size_t)bytepix;
  for (i=1; i<=naxis; i++)
  {
    snprintf((keyword+5), KEYWORD_SIZE-5, "%d", i);
    if ((p=get_FITS_key(ncards, header, keyword, val)) == -1)
    {
      printf("\n\tERROR! FITS_image_size(): %s not found in the header\n", keyword);
      return(-1);
    }
    if (p != i+2)
      printf("WARNING! header does not conform to FITS standard (%s)\n",
        keyword);

    sscanf(val, "%d", &naxisn);
    size*=naxisn;
  }

  return(size);
}
/*--------------------------------------------------------*/
/*  Read FITS header (primary or extension)               */
/*  Return number of card images.                         */
/*--------------------------------------------------------*/
int read_FITS_header(FILE *inf, char ***header)
{
        char  **lheader;
        short eoh;            // end of header
        int   i,
              hc,             // number of header cards
              hr;             // number of header records

  eoh=0;

  if ((lheader=(char **)calloc(RECORD_CARDS, sizeof(char *))) == NULL)
    errmess("read_FITS_header(): calloc(lheader)");

  hc=0;
  for (hr=1; eoh==0; hr++)
  {
    if (!(lheader=(char **)realloc(lheader, hr*RECORD_CARDS*sizeof(char *))))
      errmess("read_FITS_header(): realloc(lheader)");

    for (i=0; i<RECORD_CARDS; i++)
    {
      if (!(lheader[hc+i]=(char *)calloc(CARD_SIZE, sizeof(char))))
        errmess("read_FITS_header(): calloc(lheader[hc+i])");

      if (fread(lheader[hc+i], sizeof(char), CARD_SIZE, inf) != CARD_SIZE)
      {
        printf("\n\tERROR! read_FITS_header(): header corrupted\n");
        exit(EXIT_FAILURE);
      }
    }

    for (i=0; i<RECORD_CARDS; i++)
      if (!strncmp(lheader[hc+i], "END     ", KEYWORD_SIZE)) eoh=1;

    hc+=RECORD_CARDS;
  }

  *header=lheader;

  return(hc);
}
/*--------------------------------------------------------*/
/*  Read primary and extensions headers                   */
/*  Return number of extensions (0=primary header only)   */
/*  Byte offsets for extension headers recorded.          */
/*--------------------------------------------------------*/
int read_FITS_headers(FILE *inf, int **ncards, char ****header, size_t **offsets)
{
        char    ***lheader,         // local pointer to the header
                value[VALUE_SIZE],  // card value buffer
                val[VALUE_SIZE],    // card value buffer
                tmp,
                eoF;
        int     p,
                hen,          // number of header extensions (0=primary header)
                *lncards;     // size of the header extension
        long    offset;       // FITS file offset
        size_t  imsize,       // image size
                *loffsets;    // table of FITS file offsets to the extensions

/* read primary header */
  if (!(lheader=(char ***)calloc(1, sizeof(char **))))
    errmess("read_FITS_headers(): calloc(lheader)");
  if (!(lncards=(int *)calloc(1, sizeof(int))))
    errmess("read_FITS_headers(): calloc(lncards)");
  if (!(loffsets=(size_t *)calloc(1, sizeof(size_t))))
    errmess("read_FITS_headers(): calloc(loffsets)");

  loffsets[0]=0L;
  if (!(lncards[0]=read_FITS_header(inf, &lheader[0]))) return(-1);

/* check conformance with the FITS standard */
  if ((p=get_FITS_key(lncards[0], lheader[0], "SIMPLE", value)) == -1)
  {
    printf("\n\tERROR! read_FITS_headers(): SIMPLE not found in the header\n");
    return(-1);
  }
  if (p != 0)
    printf("WARNING! Header does not conform the FITS standard (SIMPLE)\n");
  sscanf(value, "%s", val);
  if ((strcmp(val, "T")) && (strncmp(val, "T/", 2)))
    printf("WARNING! File does not conform the FITS standard (SIMPLE)\n");

/* check for possible extensions */
  if ((p=get_FITS_key(lncards[0], lheader[0], "EXTEND", value)) == -1)
  {
    *ncards=lncards;
    *header=lheader;
    *offsets=loffsets;

    return(0);
  }

  if (p < 3)
    printf("WARNING! Header does not conform the FITS standard (EXTEND)\n");
  sscanf(value, "%s", val);
  if ((strcmp(val, "T")) && (strncmp(val, "T/", 2)))
  {
    *ncards=lncards;
    *header=lheader;
    *offsets=loffsets;

    return(0);
  }

/* read extension headers */
  hen=0;
  eoF=0;
  while (!eoF)
  {
/* skip the image */
    imsize=FITS_image_size(lncards[hen], lheader[hen]);
    offset=(long)imsize;
    if (imsize%RECORD_SIZE) offset+=(RECORD_SIZE-imsize%RECORD_SIZE);
    fseek(inf, offset, SEEK_CUR);

/* check for the existence of another extension */
    imsize=fread(&tmp, 1, 1, inf);
    if (feof(inf)) break;
    fseek(inf, -1, SEEK_CUR);

/* read another header */
    hen++;

    if (!(lheader=(char ***)realloc(lheader, (hen+1)*sizeof(char **))))
      errmess("read_FITS_headers(): realloc(lheader)");
    if (!(lncards=(int *)realloc(lncards, (hen+1)*sizeof(int))))
      errmess("read_FITS_headers(): realloc(lncards)");
    if (!(loffsets=(size_t *)realloc(loffsets, (hen+1)*sizeof(size_t))))
      errmess("read_FITS_headers(): calloc(loffsets)");

    loffsets[hen]=(size_t)ftell(inf);
    if (!(lncards[hen]=read_FITS_header(inf, &lheader[hen]))) return(-1);

/* check the XTENSION value - only IMAGE supported */
    if ((p=get_FITS_key(lncards[hen], lheader[hen], "XTENSION", value) == -1))
    {
      printf("\n\tERROR! read_FITS_headers(): XTENSION not found in the header\n");
      return(-1);
    }
    if (p != 0)
      printf("WARNING! Header does not conform the FITS standard (XTENSION)\n");

    if (strncmp(value, "'IMAGE   '", 10))
      printf("WARNING! Unsupported extension type: %s\n", value);

    if (!lncards[hen]) eoF=1;
  }

  *ncards=lncards;
  *header=lheader;
  *offsets=loffsets;

  return(hen);
}
/*--------------------------------------------------------*/
int del_header_card(int ncards, char ***header, char *keyword)
{
        char  **lheader,
              val[VALUE_SIZE];
        int   p,
              i;

  lheader=*header;

  if ((p=get_FITS_key(ncards, lheader, keyword, val)) == -1) return(ncards);

  for (i=p; i<ncards-1; i++)
    memcpy(lheader[i], lheader[i+1], CARD_SIZE);

  if ((p=get_FITS_key(ncards, lheader, "END", val)) == -1)
  {
    printf("\n\tERROR! del_header_card: header corrupted (END)\n");
    return(-1);
  }

  if (!(p%RECORD_CARDS))
  {
    ncards-=RECORD_CARDS;

    for (i=0; i<RECORD_CARDS; i++) free(lheader[ncards+i]);
    if (!(lheader=(char **)realloc(lheader, (ncards)*sizeof(char *))))
      errmess("del_header_card(): realloc(lheader)");

    *header=lheader;
  }

  return(ncards);
}
/*--------------------------------------------------------*/
void write_FITS_header(FILE *outf, int ncards, char **header)
{
        char  tmp[CARD_SIZE];
        int   i;

  for (i=0; i<ncards; i++)
  {
    if (fwrite(header[i], 1, CARD_SIZE, outf) != CARD_SIZE)
    {
      printf("\n\tERROR! writing FITS header\n");
      exit(EXIT_FAILURE);
    }
  }

  if (ncards%RECORD_CARDS)
  {
    memset(tmp, ' ', CARD_SIZE);
    for (i=0; i<RECORD_CARDS-ncards%RECORD_CARDS; i++)
    {
      if (fwrite(tmp, 1, CARD_SIZE, outf) != CARD_SIZE)
      {
        printf("\n\tERROR! writing FITS header\n");
        exit(EXIT_FAILURE);
      }
    }
  }

  return;
}
/*--------------------------------------------------------*/
size_t read_FITS_image(FILE *inf, char **header, int ncards, void **buf)
{
        void    *ptr;           /* data table pointer           */
        char    keyword[KEYWORD_SIZE+1],  /* header entry keyword */
                value[VALUE_SIZE];
        int     i,              /* loop numerator               */
                naxes,          /* image dimension              */
                *naxis,         /* table of sizes of axes       */
                bitpix,         /* bits per pixel               */
                bytepix;        /* bytes per pixel              */
        size_t  pixnum;         /* number of image pixels       */

  if (get_FITS_key(ncards, header, "NAXIS", value) == -1) return(0);
  sscanf(value, "%d", &naxes);

  if (!(naxis=(int *)calloc((size_t)naxes, sizeof(int))))
    errmess("read_FITS_image(): calloc(naxis)");

  for (i=0; i<naxes; i++)
  {
    snprintf(keyword, KEYWORD_SIZE, "NAXIS%d", i+1);
    if (get_FITS_key(ncards, header, keyword, value) == -1) return(0);
    sscanf(value, "%d", &naxis[i]);
  }

  pixnum=1;
  for (i=0; i<naxes; i++) pixnum*=naxis[i];

  free(naxis);

  if (get_FITS_key(ncards, header, "BITPIX", value) == -1) return(0);
  sscanf(value, "%d", &bitpix);
  bytepix=abs(bitpix)/8;

  if (!(ptr=(void *)malloc(pixnum*bytepix)))
    errmess("read_FITS_image(): malloc(ptr)");

  if (fread(ptr, (size_t)bytepix, pixnum, inf) != pixnum)
  {
    printf("\n\tERROR! reading FITS image\n");
    exit(EXIT_FAILURE);
  }

  if (needswap()) swap_bytes(ptr, pixnum, bytepix);

  *buf=ptr;

  return(pixnum);
}
/*--------------------------------------------------------*/
float *read_FITS_1Dfile(char *iname, char datatype,
                        int *ncards, char ***header, int *n1)
{
        void    *buf;
        char    value[VALUE_SIZE],
		**lhead;
        int     i,
		hs,
                bitpix,
                naxes,
                nx;
        size_t  nbytes;
        float   *ldata;
        double  FITS_bzero,
                FITS_bscale;
        FILE    *inf;

  if (!(inf=fopen(iname, "r")))  errmess(iname);

  hs=read_FITS_header(inf, &lhead);
  nbytes=read_FITS_image(inf, lhead, hs, &buf);

  fclose(inf);

  if (get_FITS_key(hs, lhead, "NAXIS", value) == -1)
  {
    printf("\n\tERROR! %s: NAXIS not found in the FITS header\n", iname);
    exit(EXIT_FAILURE);
  }
  sscanf(value, "%d", &naxes);
  if (naxes != 1)
  {
    printf("\n\tERROR! NAXIS= %d: %s is not 1 dimensional\n", naxes, iname);
    exit(EXIT_FAILURE);
  }

  if (get_FITS_key(hs, lhead, "NAXIS1", value) == -1)
  {
    printf("\n\tERROR! %s: NAXIS1 not found in the FITS header\n", iname);
    exit(EXIT_FAILURE);
  }
  sscanf(value, "%d", &nx);

  if (get_FITS_key(hs, lhead, "BITPIX", value) == -1)
  {
    printf("\n\tERROR! %s: BITPIX not found in the FITS header\n", iname);
    exit(EXIT_FAILURE);
  }
  sscanf(value, "%d", &bitpix);

  if (!(ldata=(float *)calloc((size_t)nx, sizeof(long))))
    errmess("read_FITS_1Dfile(): calloc(ldata)");

  if (get_FITS_key(hs, lhead, "BZERO", value) == -1)
    FITS_bzero=0.0;
  else
    sscanf(value, "%lg", &FITS_bzero);

  if (get_FITS_key(hs, lhead, "BSCALE", value) == -1)
    FITS_bscale=1.0;
  else
    sscanf(value, "%lg", &FITS_bscale);

  switch(bitpix)
  {
    case  16: if (datatype == 'u')
                for (i=0; i<nx; i++)
                  ldata[i]=(float)((unsigned short *)buf)[i];
              else
                for (i=0; i<nx; i++)
                  ldata[i]=(float)((short *)buf)[i];
              break;
    case  32: if (datatype == 'u')
                for (i=0; i<nx; i++)
                  ldata[i]=(float)((unsigned *)buf)[i];
              else
                for (i=0; i<nx; i++)
                  ldata[i]=(float)((int *)buf)[i];
              break;
    case  64: if (datatype == 'u')
                for (i=0; i<nx; i++)
                  ldata[i]=(float)((unsigned long *)buf)[i];
              else
                for (i=0; i<nx; i++)
                  ldata[i]=(float)((long *)buf)[i];
              break;
    case -32: for (i=0; i<nx; i++)
                ldata[i]=((float *)buf)[i];
              break;
    case -64: for (i=0; i<nx; i++)
                ldata[i]=(float)((double *)buf)[i];
              break;
    default:  printf("\n\tERROR! %s: bitpix= %d not supported\n",
                      iname, bitpix);
              exit(EXIT_FAILURE);
  }

  free(buf);

  if (FITS_bscale != 1.0) for (i=0; i<nx; i++) ldata[i]*=FITS_bscale;
  if (FITS_bzero != 0.0)  for (i=0; i<nx; i++) ldata[i]+=FITS_bzero;

  *n1=nx;
  *ncards=hs;
  *header=lhead;

  return(ldata);
}
/*--------------------------------------------------------*/
float *read_FITS_2D1file(char *iname, char datatype,
                        int *ncards, char ***header, int *n1, int *n2)
{
        void    *buf;
        char    value[VALUE_SIZE],
		**lhead;
        int     hs,
                bitpix,
                naxes,
                nx,
                ny;
        size_t  npix,
                i;
        float   *ldata;
        double  FITS_bzero,
                FITS_bscale;
        FILE    *inf;

  if (!(inf=fopen(iname, "r")))  errmess(iname);

  hs=read_FITS_header(inf, &lhead);
  npix=read_FITS_image(inf, lhead, hs, &buf);

  fclose(inf);

  if (get_FITS_key(hs, lhead, "NAXIS", value) == -1)
    errmess("NAXIS not found in header");
  sscanf(value, "%d", &naxes);
  if (naxes != 2) errmess("Image must be 2 dimensional");

  if (get_FITS_key(hs, lhead, "NAXIS1", value) == -1)
    errmess("NAXIS1 not found in header");
  sscanf(value, "%d", &nx);

  if (get_FITS_key(hs, lhead, "NAXIS2", value) == -1)
    errmess("NAXIS2 not found in header");
  sscanf(value, "%d", &ny);

  if (get_FITS_key(hs, lhead, "BITPIX", value) == -1)
    errmess("BITPIX not found in header");
  sscanf(value, "%d", &bitpix);

  if (get_FITS_key(hs, lhead, "BZERO", value) == -1)
    FITS_bzero=0.0;
  else
    sscanf(value, "%lg", &FITS_bzero);

  if (get_FITS_key(hs, lhead, "BSCALE", value) == -1)
    FITS_bscale=1.0;
  else
    sscanf(value, "%lg", &FITS_bscale);

/*  bytepix=abs(bitpix)/8; */

  if (!(ldata=(float *)calloc(npix, sizeof(float))))
    errmess("calloc(ldata)");

  switch(bitpix)
  {
    case  16: if (datatype == 'u')
                for (i=0; i<npix; i++)
                  ldata[i]=(float)((unsigned short *)buf)[i];
              else
                for (i=0; i<npix; i++)
                  ldata[i]=(float)((short *)buf)[i];
              break;
    case  32: if (datatype == 'u')
                for (i=0; i<npix; i++)
                  ldata[i]=(float)((unsigned *)buf)[i];
              else
                for (i=0; i<npix; i++)
                  ldata[i]=(float)((int *)buf)[i];
              break;
    case  64: if (datatype == 'u')
                for (i=0; i<npix; i++)
                  ldata[i]=(float)((unsigned long *)buf)[i];
              else
                for (i=0; i<npix; i++)
                  ldata[i]=(float)((long *)buf)[i];
              break;
    case -32: for (i=0; i<npix; i++)
                ldata[i]=((float *)buf)[i];
              break;
    case -64: for (i=0; i<npix; i++)
                ldata[i]=(float)((double *)buf)[i];
              break;
    default:  printf("bitpix=%d\n", bitpix);
              errmess("Not supported");
  }

  free(buf);

  if (FITS_bscale != 1.0)
    for (i=0; i<npix; i++)
      ldata[i]*=FITS_bscale;

  if (FITS_bzero != 0.0)
    for (i=0; i<npix; i++)
      ldata[i]+=FITS_bzero;

  *n1=nx;
  *n2=ny;
  *header=lhead;
  *ncards=hs;

  return(ldata);
}
/*--------------------------------------------------------*/
void scale(int size, float *data, short *idata)
{
        int     i;
        float   min,
                max;
        double  FITS_bscale,
                FITS_bzero;

/** scalling factors **/
  min=max=data[0];
  for (i=1; i<size; i++)
  {
    if (data[i] > max) max=data[i];
    if (data[i] < min) min=data[i];
  }
  FITS_bzero=(max-min)/2.0;
  FITS_bscale=(max-min)/65535.0;

  for (i=0; i<size; i++)
    idata[i]=(short)(nearbyint(data[i]-FITS_bzero)/FITS_bscale);

  return;
}
/*--------------------------------------------------------*/
void write_FITS_1Dimage(FILE *outf, size_t npix, int bytepix, void *data)
{
        char    *tmp;
        size_t  ls,
                size;

  if (needswap()) swap_bytes(data, npix, bytepix);
  if (fwrite(data, (size_t)bytepix, npix, outf) != npix)
  {
    printf("\n\tERROR! writing FITS image\n");
    exit(EXIT_FAILURE);
  }

  size=npix*bytepix;
  if ((ls=size%RECORD_SIZE) != 0) ls=RECORD_SIZE-ls;
  if (!(tmp=(char *)malloc(ls)))  errmess("write_FITS_1Dimage: malloc(tmp)");
  memset(tmp, 0, ls);
  if (fwrite(tmp, 1, ls, outf) != ls)
  {
    printf("\n\tERROR! writing FITS image\n");
    exit(EXIT_FAILURE);
  }

  free(tmp);

  return;
}
/*--------------------------------------------------------*/
void write_FITS_1Dfile(char *oname, int ncards, char **header,
                        size_t npix, int bytepix, void *data)
{
        FILE  *outf;

  if (!(outf=fopen(oname, "w"))) errmess(oname);

  write_FITS_header(outf, ncards, header);
  write_FITS_1Dimage(outf, npix, bytepix, data);

  fclose(outf);

  return;
}
/*--------------------------------------------------------*/
void write_FITS_2D1file(char *oname, int ncards, char **header,
                        size_t npix, int bytepix, void *data)
{
  write_FITS_1Dfile(oname, ncards, header, npix, bytepix, data);

  return;
}
/*** END ***/
