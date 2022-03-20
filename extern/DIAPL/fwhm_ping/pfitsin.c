/*========================================================*/
/*                                                        */
/*  pfitsin.c       version 3.1.3         2006.02.14      */
/*                                                        */
/*  Copyright (C) 2006 by Wojtek Pych, CAMK PAN           */
/*                                                        */
/*  Written for GNU project C and C++ Compiler            */
/*  Private FITS I/O library.                             */
/*                                                        */
/*--------------------------------------------------------*/
/*                                                        */
/* Subset of pfitsio.c functions.                         */
/*                                                        */
/*========================================================*/

#include "pfitsin.h"
#include "errmess.h"

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
void swap_bytes(void *ptr, int ndata, int nbytes)
{
        char  c,
              *data;
        int   i, j, k, m;

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
      errmess("calloc(loffsets)");

    loffsets[hen]=offset;
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
int read_FITS_image(FILE *inf, char **header, int hsize, void **buf)
{
        void  *ptr;           /* data table pointer           */
        char  flag[9],        /* header entry flag            */
              hline[72];
        int   i,              /* loop numerator               */
              naxes,          /* image dimension              */
              *naxis,         /* table of sizes of axes       */
              pixnum,         /* number of image pixels       */
              bitpix,
              bytepix;

  if (get_FITS_key(hsize, header, "NAXIS", hline) == -1) return(0);
  sscanf(hline, "%d", &naxes);

  if (!(naxis=(int *)calloc(naxes, sizeof(int))))
    errmess("pfitsio: read_FITS_image - calloc(naxis)");

  for (i=0; i<naxes; i++)
  {
    sprintf(flag, "NAXIS%d", i+1);
    if (get_FITS_key(hsize, header, flag, hline) == -1) return(0);
    sscanf(hline, "%d", &naxis[i]);
  }

  pixnum=1;
  for (i=0; i<naxes; i++) pixnum*=naxis[i];

  free(naxis);

  if (get_FITS_key(hsize, header, "BITPIX", hline) == -1) return(0);
  sscanf(hline, "%d", &bitpix);
  bytepix=abs(bitpix)/8;

  if (!(ptr=(void *)malloc(pixnum*bytepix)))
    errmess("pfitsio: read_FITS_image - malloc(ptr)");

  if (fread(ptr, bytepix, pixnum, inf) != pixnum)
    errmess("pfitsio: read_FITS_image - reading data");

  if (needswap()) swap_bytes(ptr, pixnum, bytepix);

  *buf=ptr;

  return(pixnum);
}
/*--------------------------------------------------------*/
float *read_FITS_1Dfile(char *iname, char datatype, int *hsize, char ***header,
                        int *n1)
{
        void    *buf;
        char    hline[72],
		**lhead;
        int     i,
		hs,
                bitpix,
                naxes,
                nx;
        float   *ldata;
        double  FITS_bzero,
                FITS_bscale;
        FILE    *inf;

  if (!(inf=fopen(iname, "r")))  errmess(iname);

  hs=read_FITS_header(inf, &lhead);
  read_FITS_image(inf, lhead, hs, &buf);

  fclose(inf);

  if (get_FITS_key(hs, lhead, "NAXIS", hline) == -1)
    errmess("pfitsio: NAXIS not found in header");
  sscanf(hline, "%d", &naxes);
  if (naxes != 1) errmess("Image must be 1 dimensional");

  if (get_FITS_key(hs, lhead, "NAXIS1", hline) == -1)
    errmess("pfitsio: NAXIS1 not found in header");
  sscanf(hline, "%d", &nx);

  if (get_FITS_key(hs, lhead, "BITPIX", hline) == -1)
    errmess("pfitsio: BITPIX not found in header");
  sscanf(hline, "%d", &bitpix);

  if (!(ldata=(float *)calloc(nx, sizeof(long))))
    errmess("calloc(ldata)");

  if (get_FITS_key(hs, lhead, "BZERO", hline) == -1)
    FITS_bzero=0.0;
  else
    sscanf(hline, "%lg", &FITS_bzero);

  if (get_FITS_key(hs, lhead, "BSCALE", hline) == -1)
    FITS_bscale=1.0;
  else
    sscanf(hline, "%lg", &FITS_bscale);

  switch(bitpix)
  {
    case   8: if (datatype == 'u')
                for (i=0; i<nx; i++)
                  ldata[i]=(float)((unsigned char *)buf)[i];
              else
                for (i=0; i<nx; i++)
                  ldata[i]=(float)((char *)buf)[i];
              break;
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
    case -32: for (i=0; i<nx; i++)
                ldata[i]=((float *)buf)[i];
              break;
    case -64: for (i=0; i<nx; i++)
                ldata[i]=(float)((double *)buf)[i];
              break;
    default:  printf("bitpix=%d\n", bitpix);
              errmess("Not supported");
  }

  free(buf);

  if (FITS_bscale != 1.0) for (i=0; i<nx; i++) ldata[i]*=FITS_bscale;
  if (FITS_bzero != 0.0)  for (i=0; i<nx; i++) ldata[i]+=FITS_bzero;

  *n1=nx;
  *header=lhead;
  *hsize=hs;

  return(ldata);
}
/*--------------------------------------------------------*/
float **read_FITS_2Dfile(char *iname, char datatype,
                         int *hsize, char ***header, int *n1, int *n2)
{
        void    *buf;
        char    hline[72],
		**lhead;
        int     i,
                j,
		hs,
                bitpix,
                naxes,
                nx,
                ny;
        float   **ldata;
        double  FITS_bzero,
                FITS_bscale;
        FILE    *inf;

  if (!(inf=fopen(iname, "r")))  errmess(iname);

  hs=read_FITS_header(inf, &lhead);
  read_FITS_image(inf, lhead, hs, &buf);

  fclose(inf);

  if (get_FITS_key(hs, lhead, "NAXIS", hline) == -1)
    errmess("pfitsio: NAXIS not found in header");
  sscanf(hline, "%d", &naxes);
  if (naxes != 2) errmess("Image must be 2 dimensional");

  if (get_FITS_key(hs, lhead, "NAXIS1", hline) == -1)
    errmess("pfitsio: NAXIS1 not found in header");
  sscanf(hline, "%d", &nx);

  if (get_FITS_key(hs, lhead, "NAXIS2", hline) == -1)
    errmess("pfitsio: NAXIS2 not found in header");
  sscanf(hline, "%d", &ny);

  if (get_FITS_key(hs, lhead, "BITPIX", hline) == -1)
    errmess("pfitsio: BITPIX not found in header");
  sscanf(hline, "%d", &bitpix);

/*  bytepix=abs(bitpix)/8; */

  if (!(ldata=(float **)calloc(ny, sizeof(float *))))
    errmess("calloc(ldata)");
  for (i=0; i<ny; i++)
  {
    if (!(ldata[i]=(float *)calloc(nx, sizeof(float))))
      errmess("calloc(ldata[i])");
  }

  if (get_FITS_key(hs, lhead, "BZERO", hline) == -1)
    FITS_bzero=0.0;
  else
    sscanf(hline, "%lg", &FITS_bzero);

  if (get_FITS_key(hs, lhead, "BSCALE", hline) == -1)
    FITS_bscale=1.0;
  else
    sscanf(hline, "%lg", &FITS_bscale);

  switch(bitpix)
  {
    case   8: if (datatype == 'u')
                for (i=0; i<ny; i++)
                  for (j=0; j<nx; j++)
                    ldata[i][j]=(float)((unsigned char *)buf)[i*nx+j];
              else
                for (i=0; i<ny; i++)
                  for (j=0; j<nx; j++)
                    ldata[i][j]=(float)((char *)buf)[i*nx+j];
              break;
    case  16: if (datatype == 'u')
                for (i=0; i<ny; i++)
                  for (j=0; j<nx; j++)
                    ldata[i][j]=(float)((unsigned short *)buf)[i*nx+j];
              else
                for (i=0; i<ny; i++)
                  for (j=0; j<nx; j++)
                    ldata[i][j]=(float)((short *)buf)[i*nx+j];
              break;
    case  32: if (datatype == 'u')
                for (i=0; i<ny; i++)
                  for (j=0; j<nx; j++)
                    ldata[i][j]=(float)((unsigned *)buf)[i*nx+j];
              else
                for (i=0; i<ny; i++)
                  for (j=0; j<nx; j++)
                    ldata[i][j]=(float)((int *)buf)[i*nx+j];
              break;
    case -32: for (i=0; i<ny; i++)
                for (j=0; j<nx; j++)
                  ldata[i][j]=((float *)buf)[i*nx+j];
              break;
    case -64: for (i=0; i<ny; i++)
                for (j=0; j<nx; j++)
                  ldata[i][j]=(float)((double *)buf)[i*nx+j];
              break;
    default:  printf("bitpix=%d\n", bitpix);
              errmess("Not supported");
  }

  free(buf);

  if (FITS_bscale != 1.0)
    for (i=0; i<ny; i++)
      for (j=0; j<nx; j++)
        ldata[i][j]*=FITS_bscale;

  if (FITS_bzero != 0.0)
    for (i=0; i<ny; i++)
      for (j=0; j<nx; j++)
        ldata[i][j]+=FITS_bzero;

  *n1=nx;
  *n2=ny;
  *header=lhead;
  *hsize=hs;

  return(ldata);
}
/*** END ***/
