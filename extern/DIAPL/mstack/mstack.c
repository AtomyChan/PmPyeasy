/*========================================================*/
/*                                                        */
/*  mstack.c            version 1.9.3   2006.10.24        */
/*                                                        */
/*  Original source code:                                 */
/*  Copyright (C) 2006 by Przemek Wozniak                 */
/*  wozniak@lanl.gov                                      */
/*                                                        */
/*  Modifications:                                        */
/*  Copyright (C) 2006 by Wojtek Pych, CAMK PAN           */
/*  pych@camk.edu.pl                                      */
/*                                                        */
/*  Written for GNU project C and C++ Compiler            */
/*                                                        */
/* Program stacks images using median for a given pixel   */
/* and saving pixels which were bad on some images,       */
/* but have enough information in them.                   */
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

#include <stdarg.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "defs.h"
#include "funcs.h"
#include "errmess.h"
#include "pfitsio1.h"

#define BIG_FLOAT 1.0e6
#define NAME_LEN 200

/*--------------------------------------------------------*/
void usage(void)
{
  printf("\n\tUSAGE: mstack parameter_file instrument_file ");
  printf(" reference_image image_names_file output_image\n\n");
  exit(EXIT_FAILURE);
}
/*--------------------------------------------------------*/
int read_imnames(char *iname, char ***imnames)
{
        char  **names,
              line[NAME_LEN];
        int   i;
        FILE  *inpf;

  if (!(names=(char **)calloc(1, sizeof(char *))))
    errmess("read_imnames: calloc(names)");

  if (!(inpf=fopen(iname, "r"))) errmess(iname);

  for (i=0; !feof(inpf); i++)
  {
    if (fgets(line, NAME_LEN, inpf) == NULL) break;
    if (feof(inpf)) break;

    if (!(names=(char **)realloc(names, (i+1)*sizeof(char *))))
      errmess("read_imnames: realloc(names)");
    if (!(names[i]=(char *)calloc(strlen(line)+1, sizeof(char))))
      errmess("read_imnames: calloc(names[i])");

    sscanf(line, "%s", names[i]);
  }

  fclose(inpf);

  *imnames=names;

  return(i);
}
/*--------------------------------------------------------*/
void hedit(int ncards, char **header, char *name, char type, char *descr,
           char *fmt, ...)
{
        char    cval, *sval, entry[VALUE_SIZE];
        int     len, i, ival, found;
        double  dval;
        va_list ap;

  va_start(ap, fmt);

  switch(type)
  {
    case 'c': cval = (char)va_arg(ap, int);
              len  = snprintf(entry, VALUE_SIZE, fmt, cval);
              break;
    case 's': sval = va_arg(ap, char *);
              len  = snprintf(entry, VALUE_SIZE, fmt, sval);
              break;
    case 'i': ival = va_arg(ap, int);
              len  = snprintf(entry, VALUE_SIZE, fmt, ival);
              break;
    case 'd':
    case 'f': dval = va_arg(ap, double);
              len  = snprintf(entry, VALUE_SIZE, fmt, dval);
              break;
    default:  printf("hedit: warning: unknown type of keyword - ");
              printf("header unchanged\n");
              len = 0;
              break;
  }

  if (len > 0)
  {
    if (len < VALUE_SIZE-6)
      len += snprintf(entry+len, (size_t)(VALUE_SIZE-len), " / %s", descr);

    for (i=len; i<VALUE_SIZE; i++) entry[i] = ' ';

    found=0;
    for (i=0; i<ncards; i++)
    {
      if (strncmp(header[i], name, KEYWORD_SIZE) == 0)
      {
        memcpy(header[i]+10, entry, VALUE_SIZE);
        found = 1;
      }
    }

    if (!found)
      printf("hedit: warning: %s not found - value unchanged\n", name);
  }

  return;
}
/*--------------------------------------------------------*/
void sort(int n, float *data)
{
        int   i,
              j;
        float tmp;

  for (i=0; i<n-1; i++)
  {
    for (j=i+1; j<n; j++)
    {
      if (data[i] > data[j])
      {
        tmp=data[i];
        data[i]=data[j];
        data[j]=tmp;
      }
    }
  }

  return;
}
/*--------------------------------------------------------*/
float fcalc_median(int size, float *data)
{
        float median;

  if (size%2) median=data[size/2];
  else        median=(data[size/2]+data[size/2-1])/2.0;

  return(median);
}
/*--------------------------------------------------------*/
float specstat(int n1, int n2, float *data, float *stddev)
{
        int     i,
                n;
        double  mean, lstddev;

  mean = lstddev = 0.0;
  n=n2-n1+1;
  for (i=n1; i<n2; i++)
  {
    mean+=data[i]/n;
    lstddev+=data[i]/n*data[i];
  }

  *stddev=(float)sqrt(lstddev-mean*mean);

  return((float)mean);
}
/*--------------------------------------------------------*/
void  stack_pix(int npix, int nim, float **data, float *bkg, float *scale,
                float refbg, float *imstack, PAR_STRUCT par)
{
        int   i,    // loop over images
              j,    // loop over pixels
              k;    // number of good pixels
        float *stat,
              mean,
              stddev;

  if (!(stat=(float *)calloc((size_t)nim, sizeof(float))))
    errmess("stack_pix: calloc(stat)");

  for (j=0; j<npix; j++)
  {
    k=0;

    for (i=0; i<nim; i++)
    {
      if (data[i][j] > par.sat_level)
      {
        if (par.verbose > 5)
          printf("data[%d](%d, %d)= %g > %g (SAT_LEVEL)\n",
                  i, j%par.nx0+1, j/par.nx0+1, data[i][j], par.sat_level);
        continue;
      }

      if (data[i][j] < par.min_level)
      {
        if (par.verbose > 5)
          printf("data[%d](%d, %d)= %g < %g (MIN_LEVEL)\n",
                  i, j%par.nx0+1, j/par.nx0+1, data[i][j], par.min_level);
        continue;
      }

      stat[k++]=(data[i][j]-bkg[i])*scale[i];
    }

    if (k < par.min_ngood)
    {
      if (par.verbose > 4)
        printf("stack_pix:  (%d, %d) num. of valid data= %d < %d (MIN_NGOOD)\n",
                j%par.nx0+1, j/par.nx0+1, k, par.min_ngood);
      imstack[j]=par.bad_value;
      continue;
    }

    sort(k, stat);

    if (par.median)
    {
      imstack[j]=fcalc_median(k, stat)+refbg;
      continue;
    }

    mean=specstat(par.nlow, k-par.nhigh, stat, &stddev);

    if (par.verbose > 5)
      printf("k= %d  mean= %g +- %g\n", k-par.nlow-par.nhigh, mean, stddev);

    imstack[j]=mean+refbg;
  }

  free(stat);

  return;
}
/*--------------------------------------------------------*/
int modify_header(int ncards, char ***header)
{
        char  **lheader,
              value[VALUE_SIZE],
              card[CARD_SIZE];
        int   k;

  lheader=*header;

  if ((k=get_FITS_key(ncards, lheader, "BITPIX", value)) == -1)
  {
    printf("\n\tERROR! mstack: BITPIX not found in the FITS header\n");
    exit(EXIT_FAILURE);
  }

  memset(card, ' ', CARD_SIZE);
  memcpy(card, "BITPIX  =                  -32 / Bits per pixel", 47);
  memcpy(lheader[k], card, CARD_SIZE);

  ncards=del_header_card(ncards, &lheader, "BZERO");
  ncards=del_header_card(ncards, &lheader, "BSCALE");

  *header=lheader;

  return(ncards);
}
/*--------------------------------------------------------*/
int main(int argc, char *argv[])
{
        char        **header,
                    **theader,
                    value[VALUE_SIZE],
                    **imnames,  /* input files names              */
                    *outfname,  /* output file name               */
                    *parfname,  /* parameters file name pointer   */
                    *instrname,
                    *imnfname,  /* list of input files file name  */
                    *refname;   /* background and scale reference image name */
        int         i,
                    nx0, ny0,
                    nx, ny,
                    npix,
                    k, iim, nim,
                    x0_off, y0_off,
                    isec, nsec, igood,
                    ncards, tncards;
        long        *headlen;
        size_t      outheadlen;
        float       *hist1, *hist2,
                    *bins1, *bins2,
                    low, high, low0, high0,
                    bin_size,
                    *bkg, *scale, bg, refbg, sum_xy, sum_x2,
                    *im, *refim, *imstack, **data;
        PAR_STRUCT  par;
        FILE        *outf;

/* IO stuff */
  if (argc != 6) usage();

  parfname =argv[1];
  instrname=argv[2];
  refname  =argv[3];
  imnfname =argv[4];
  outfname =argv[5];

/* read in parameters */
  get_params(parfname, instrname, &par);

  if (par.verbose > 1)
  {
    printf("parfname= %s\n", parfname);
    printf("instrname= %s\n", instrname);
    printf("refname= %s\n", refname);
    printf("imnfname= %s\n", imnfname);
    printf("outfname= %s\n", outfname);
  }

  nim=read_imnames(imnfname, &imnames);

  if (par.verbose)
  {
    printf("%s: %d images\n", imnfname, nim);
    if (par.verbose > 1)
      for (i=0; i<nim; i++) printf("%s\n", imnames[i]);
    printf("\n");
  }

  if (nim < par.min_ngood)
  {
    printf("\n\tERROR! MIN_NGOOD= %d larger than the number of images (%d)\n",
            par.min_ngood, nim);
    exit(EXIT_FAILURE);
  }

  if (nim <= par.nlow+par.nhigh)
  {
    printf("\n\tERROR! N_REJ_LOW+N_REJ_HIGH larger than the number of images\n");
    exit(EXIT_FAILURE);
  }

/* make histogram axes */
  if (!(hist1=(float *)malloc(par.nbin*sizeof(float))))
    errmess("malloc(hist1)");
  if (!(hist2=(float *)malloc((par.nbin/par.nsmooth)*sizeof(float))))
    errmess("malloc(hist2)");
  if (!(bins1=(float *)malloc(par.nbin*sizeof(float))))
    errmess("malloc(bins1)");
  if (!(bins2=(float *)malloc(par.nbin*sizeof(float))))
    errmess("malloc(bins2)");

  bin_size = (par.high - par.low)/(float)par.nbin;
  if (par.verbose > 2) printf("histogram bin_size= %g\n", bin_size);

  for (k=0; k<par.nbin; k++)
    bins1[k] = par.low + (k+0.5)*bin_size;

  for (k=0; k<par.nbin/par.nsmooth; k++)
    bins2[k] = par.low + (k+0.5)*par.nsmooth*bin_size;

  if (par.verbose > 2) printf("bkg and scale reference: %s\n", refname);

  refim=read_FITS_2D1file(refname, 's', &ncards, &header, &nx0, &ny0);

  if (get_FITS_key(ncards, header, "BITPIX", value) == -1)
  {
    printf("\n\tERROR! mstack: BITPIX not found in the header (%s)\n", refname);
    exit(EXIT_FAILURE);
  }
  if (atoi(value) != -32)
  {
    printf("\n\tERROR! %s: image should have BITPIX = -32\n", refname);
    exit(EXIT_FAILURE);
  }
  if (par.verbose > 3) printf("%s: %d x %d\n", refname, nx0, ny0);

  if (nx0 != par.nx0)
  {
    printf("\n\tERROR! %s: NAXIS1= %d <> NX0= %d\n", refname, nx0, par.nx0);
    printf("mstack: wrong image size\n");
    exit(EXIT_FAILURE);
  }
  if (ny0 != par.ny0)
  {
    printf("\n\tERROR! %s: NAXIS2= %d <> NY0= %d\n", refname, ny0, par.ny0);
    printf("mstack: wrong image size\n");
    exit(EXIT_FAILURE);
  }

  histogram(refim, nx0, ny0, hist1, par.nbin, par.low, par.high);
  if (par.verbose > 3) printf("Histogram created\n");

  bin(hist1, par.nbin, hist2, par.nsmooth);

  get_peak(hist2, bins2, 2, par.nbin/par.nsmooth-2, &low0, &high0);

  refbg = median(hist1, bins1, par.nbin, low0, high0, par.bkg_frac);

  if (par.verbose)
    printf("reference image: %s  bkg= %g  scale= 1\n", refname, refbg);

/* get memory */
  if (!(bkg  =(float *)malloc(nim*sizeof(float)))) errmess("malloc(bkg)");
  if (!(scale=(float *)malloc(nim*sizeof(float)))) errmess("malloc(scale)");
  if (!(headlen=(long *)malloc(nim*sizeof(long)))) errmess("malloc(headlen)");

  igood = 0;

/**********************************************************/
/* loop over remaining images to get background and scale */
/**********************************************************/

  for (iim=0; iim<nim; iim++)
  {
    if (par.verbose > 3) printf("image[%d]: %s\n", iim, imnames[iim]);
    im=read_FITS_2D1file(imnames[iim], 's', &tncards, &theader, &nx0, &ny0);

    for (i=0; i<tncards; i++) free(theader[i]);
    free(theader);

    headlen[iim]=tncards*CARD_SIZE;
    if (tncards%RECORD_CARDS)
      headlen[iim]+=(RECORD_CARDS-tncards%RECORD_CARDS)*CARD_SIZE;

    if (nx0 != par.nx0)
    {
      printf("\n\tERROR! %s: NAXIS1= %d <> NX0= %d\nmstack: wrong image size\n",
              imnames[iim], nx0, par.nx0);
      exit(EXIT_FAILURE);
    }

    if (ny0 != par.ny0)
    {
      printf("\n\tERROR! %s: NAXIS2= %d <> NY0= %d\nmstack: wrong image size\n",
              imnames[iim], ny0, par.ny0);
      exit(EXIT_FAILURE);
    }

    histogram(im, nx0, ny0, hist1, par.nbin, par.low, par.high);

    bin(hist1, par.nbin, hist2, par.nsmooth);

    get_peak(hist2, bins2, 2, par.nbin/par.nsmooth-2, &low, &high);

    bg=median(hist1, bins1, par. nbin, low, high, par.bkg_frac);

    sum_xy = sum_x2 = 0.0;

    for (k=0; k<nx0*ny0; k++)
    {
      if (refim[k] >= high0 + par.threshold)
      {
        sum_x2 += (im[k]-bg)*(im[k]-bg)/(im[k]+refim[k]);
        sum_xy += (im[k]-bg)*(refim[k]-refbg)/(im[k]+refim[k]);
      }
    }

    bkg[igood]   = bg;
    scale[igood] = sum_xy/sum_x2;

    if (par.verbose)
      printf("image: %s  bkg= %g  scale= %g     ",
              imnames[iim], bkg[igood], scale[igood]);

    if ((scale[igood] >= par.min_scale) && (scale[igood] <= par.max_scale))
    {
      if (igood != iim) strcpy(imnames[igood], imnames[iim]);
      igood++;
      if (par.verbose) printf("\n");
    }
    else  printf("rejected\n");

    free(im);
  }

  nim = igood;

  if (nim < 1)
  {
    printf("\n\tERROR! mstack: no images left!\n");
    exit(EXIT_FAILURE);
  }

/* set sizes and shifts for reading image sections */
  nx=nx0;
  ny=par.ny;

  x0_off=0;

  nsec=par.ny0/par.ny;

/* get memory for image sections */
  npix=nx*ny;
  if (!(imstack=(float *)malloc(npix*sizeof(float))))
    errmess("malloc(imstack)");

  if (!(data=(float **)malloc(nim*sizeof(float *)))) errmess("malloc(data)");
  for(iim=0; iim<nim; iim++)
    if (!(data[iim]=(float *)malloc(nx*ny*sizeof(float))))
      errmess("malloc(data[iim])");

/* modify header */
  ncards=modify_header(ncards, &header);

/* write a random output image */
  if (par.verbose > 3)  printf("\nWriting header to %s\n", outfname);
  if (!(outf=fopen(outfname, "w"))) errmess(outfname);
  write_FITS_header(outf, ncards, header);
  fclose(outf);

  for (i=0; i<ncards; i++) free(header[i]);
  free(header);

/***************************************/
/*** loop over image sections in "y" ***/
/***************************************/

  outheadlen=ncards*CARD_SIZE;

  for (isec=0; isec<=nsec; isec++)
  {
    y0_off = isec*ny;

/* leftover section may have a different size */
    if (isec == nsec) ny = ny0 - nsec*ny;
    if (ny <= 0) break;

/* read a single section of all images */
    for (iim=0; iim<nim; iim++)
      read_sector(imnames[iim], (size_t)headlen[iim], nx0, ny0, x0_off, y0_off,
                  (char *)data[iim], (int)sizeof(float), nx, ny);

    stack_pix(npix, nim, data, bkg, scale, refbg, imstack, par);

    write_sector(outfname, outheadlen, nx0, ny0, x0_off, y0_off,
                 (char *)imstack, (int)sizeof(float), nx, ny, 0);
  }

  if (par.verbose)  printf("\nimage stack written to %s\n", outfname);

  for (i=0; i<nim; i++) free(imnames[i]);
  free(imnames);

  free(bkg);
  free(scale);
  free(refim);
  free(hist1);
  free(hist2);
  free(bins1);
  free(bins2);
  free(headlen);

  return(0);
}
/*** END ***/
