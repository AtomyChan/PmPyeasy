/*========================================================*/
/*                                                        */
/*  pfitshead.c     version 3.3.3         2006.02.13      */
/*                                                        */
/*  Copyright (C) 2006 by Wojtek Pych, CAMK PAN           */
/*                                                        */
/*  Written for GNU project C and C++ Compiler            */
/*  Private FITS header libary.                           */
/*                                                        */
/*--------------------------------------------------------*/
/*                                                        */
/*  Latest modifications:                                 */
/*  FITS Extensions added.                                */
/*  get_FITS_key()                                        */
/*                                                        */
/*========================================================*/

#include "pfitshead.h"
#include "errmess.h"

/*--------------------------------------------------------*/
/*  Read header card value with comment                   */
/*  Return card number (starting from 0)                  */
/*  or -1 on error.                                       */
/*--------------------------------------------------------*/
int get_FITS_key(int hsize, char **header, char *keyword, char *value)
{
        char  lkeyword[KEYWORD_SIZE];
        int   i;              /*  loop numerator          */

  if (strlen(keyword) > KEYWORD_SIZE)
  {
    printf("WARNING! get_FITS_key(%s): invalid keyword\n", keyword);
    return(-1);
  }

  for (i=0; i<strlen(keyword); i++)             lkeyword[i]=keyword[i];
  for (i=strlen(keyword); i<KEYWORD_SIZE; i++)  lkeyword[i]=' ';

  for (i=0; i<hsize; i++)
  {
    if (!strncmp(header[i], lkeyword, KEYWORD_SIZE))
    {
      memcpy(value, header[i]+10, VALUE_SIZE);  // value may start in column 11
      break;
    }
  }

  if (i == hsize) return(-1);

  return(i);
}
/*--------------------------------------------------------*/
/*  Return: image size in bytes (without the filling)     */
/*          0 no image                                    */
/*         -1 error                                       */
/*--------------------------------------------------------*/
long FITS_image_size(int hsize, char **header)
{
        char  val[VALUE_SIZE],
              keyword[KEYWORD_SIZE];
        int   i,
              p,
              bitpix,
              bytepix,
              naxis,
              naxisn;
        long  size;

  if ((p=get_FITS_key(hsize, header, "BITPIX", val)) == -1)
  {
    printf("ERROR! FITS_image_size(): BITPIX not found in the header\n");
    return(-1);
  }
  if (p != 1)
    printf("WARNING! header does not conform to FITS standard (BITPIX)\n");
  sscanf(val, "%d", &bitpix);
  bytepix=abs(bitpix)/8;

  strcpy(keyword, "NAXIS");
  if ((p=get_FITS_key(hsize, header, keyword, val)) == -1)
  {
    printf("ERROR! FITS_image_size(): NAXIS not found in the header\n");
    return(-1);
  }
  if (p != 2)
    printf("WARNING! header does not conform to FITS standard (NAXIS)\n");

  sscanf(val, "%d", &naxis);
  if (naxis == 0) return(0L);

  size=bytepix;
  for (i=1; i<=naxis; i++)
  {
    sprintf((keyword+5), "%d", i);
    if ((p=get_FITS_key(hsize, header, keyword, val)) == -1)
    {
      printf("ERROR! FITS_image_size(): %s not found in the header\n", keyword);
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
        char  **lheader,
              eoh;            // end of header
        int   i,
              hc,             // number of header cards
              hr;             // number of header records

  eoh=0;

  if (!(lheader=(char **)calloc(RECORD_CARDS, sizeof(char *))))
    errmess("read_FITS_header(): calloc(lheader)");

  hc=0;
  for (hr=1; !eoh; hr++)
  {
    if (!(lheader=(char **)realloc(lheader, hr*RECORD_CARDS*sizeof(char *))))
      errmess("read_FITS_header(): realloc(lheader)");

    for (i=0; i<RECORD_CARDS; i++)
    {
      if (!(lheader[hc+i]=(char *)calloc(CARD_SIZE, sizeof(char))))
        errmess("read_FITS_header(): calloc(lheader[hc+i])");

      if (fread(lheader[hc+i], sizeof(char), CARD_SIZE, inf) != CARD_SIZE)
      {
        printf("ERROR! read_FITS_header(): header corrupted\n");
        return(0);
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
int read_FITS_headers(FILE *inf, int **hsize, char ****header, long **offsets)
{
        char  ***lheader,
              value[VALUE_SIZE],
              val[VALUE_SIZE],
              tmp,
              eoF;
        int   p,
              hen,            // header extension number (0=primary header)
              *lhsize;        // size of the header extension
        long  imsize,
              offset,
              *loffsets;

/* read primary header */
  if (!(lheader=(char ***)calloc(1, sizeof(char **))))
    errmess("read_FITS_headers(): calloc(lheader)");
  if (!(lhsize=(int *)calloc(1, sizeof(int))))
    errmess("read_FITS_headers(): calloc(lhsize)");
  if (!(loffsets=(long *)calloc(1, sizeof(long))))
    errmess("read_FITS_headers(): calloc(loffsets)");

  loffsets[0]=0L;
  if (!(lhsize[0]=read_FITS_header(inf, &lheader[0]))) return(-1);

/* check conformance with the FITS standard */
  if ((p=get_FITS_key(lhsize[0], lheader[0], "SIMPLE", value)) == -1)
  {
    printf("ERROR! read_FITS_headers(): SIMPLE not found in the header\n");
    return(-1);
  }
  if (p != 0)
    printf("WARNING! Header does not conform the FITS standard (SIMPLE)\n");
  sscanf(value, "%s", val);
  if ((strcmp(val, "T")) && (strncmp(val, "T/", 2)))
    printf("WARNING! File does not conform the FITS standard (SIMPLE)\n");

/* check for possible extensions */
  if ((p=get_FITS_key(lhsize[0], lheader[0], "EXTEND", value)) == -1)
  {
    *hsize=lhsize;
    *header=lheader;
    *offsets=loffsets;

    return(0);
  }

  if (p < 3)
    printf("WARNING! Header does not conform the FITS standard (EXTEND)\n");
  sscanf(value, "%s", val);
  if ((strcmp(val, "T")) && (strncmp(val, "T/", 2)))
  {
    *hsize=lhsize;
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
    imsize=FITS_image_size(lhsize[hen], lheader[hen]);
    offset=imsize;
    if (imsize%RECORD_SIZE) offset+=(RECORD_SIZE-imsize%RECORD_SIZE);
    fseek(inf, offset, SEEK_CUR);

/* check for the existence of another extension */
    fread(&tmp, 1, 1, inf);
    if (feof(inf)) break;
    fseek(inf, -1, SEEK_CUR);

/* read another header */
    hen++;

    if (!(lheader=(char ***)realloc(lheader, (hen+1)*sizeof(char **))))
      errmess("read_FITS_headers(): realloc(lheader)");
    if (!(lhsize=(int *)realloc(lhsize, (hen+1)*sizeof(int))))
      errmess("read_FITS_headers(): realloc(lhsize)");
    if (!(loffsets=(long *)realloc(loffsets, (hen+1)*sizeof(long))))
      errmess("read_FITS_headers(): calloc(loffsets)");

    loffsets[hen]=ftell(inf);
    if (!(lhsize[hen]=read_FITS_header(inf, &lheader[hen]))) return(-1);

/* check the XTENSION value - only IMAGE supported */
    if ((p=get_FITS_key(lhsize[hen], lheader[hen], "XTENSION", value) == -1))
    {
      printf("ERROR! read_FITS_headers(): XTENSION not found in the header\n");
      return(-1);
    }
    if (p != 0)
      printf("WARNING! Header does not conform the FITS standard (XTENSION)\n");

    if (strncmp(value, "'IMAGE   '", 10))
      printf("WARNING! Unsupported extension type: %s\n", value);

    if (!lhsize[hen]) eoF=1;
  }

  *hsize=lhsize;
  *header=lheader;
  *offsets=loffsets;

  return(hen);
}
/*--------------------------------------------------------*/
int del_header_card(int hsize, char ***header, char *keyword)
{
        char  **lheader,
              val[VALUE_SIZE];
        int   p,
              i;

  lheader=*header;

  if ((p=get_FITS_key(hsize, lheader, keyword, val)) == -1) return(hsize);

  for (i=p; i<hsize-1; i++)
    memcpy(lheader[i], lheader[i+1], CARD_SIZE);

  if ((p=get_FITS_key(hsize, lheader, "END", val)) == -1)
  {
    printf("ERROR! del_header_card: header corrupted (END)\n");
    return(-1);
  }

  if (!(p%RECORD_CARDS))
  {
    hsize-=RECORD_CARDS;

    for (i=0; i<RECORD_CARDS; i++) free(lheader[hsize+i]);
    if (!(lheader=(char **)realloc(lheader, (hsize)*sizeof(char *))))
      errmess("del_header_card(): realloc(lheader)");

    *header=lheader;
  }

  return(hsize);
}
/*--------------------------------------------------------*/
void write_FITS_header(FILE *outf, int hsize, char **header)
{
        char  tmp[CARD_SIZE];
        int   i;

  for (i=0; i<hsize; i++)
    fwrite(header[i], 1, CARD_SIZE, outf);

  if (hsize%RECORD_CARDS)
  {
    memset(tmp, ' ', CARD_SIZE);
    for (i=0; i<RECORD_CARDS-hsize%RECORD_CARDS; i++)
      fwrite(tmp, 1, CARD_SIZE, outf);
  }

  return;
}
/*** END ***/
