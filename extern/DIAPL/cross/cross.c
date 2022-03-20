/*========================================================*/
/*                                                        */
/*  cross.c             version 1.4.0   2005.11.23        */
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
#include <math.h>

#define NAME_LEN 201

#include "defs.h"
#include "funcs.h"

#include "errmess.h"
#include "pfitsin1.h"

/*--------------------------------------------------------*/
void usage()
{
    printf("\n  USAGE: cross parameter_file instrument_file output_file ");
    printf("image1 (image2 or -f image_names_file)\n\n");
    exit(1);
}
/*--------------------------------------------------------*/
int readlist(char *iname, char ***imnames)
{
        char  line[NAME_LEN],
              **names;
        int   i;
        FILE  *inf;

  if (!(names=(char **)calloc(1, sizeof(char *)))) errmess("calloc(names)");

  if (!(inf=fopen(iname, "r"))) errmess(iname);

  for (i=0; !feof(inf); i++)
  {
    fgets(line, NAME_LEN, inf);
    if (feof(inf)) break;

    if (!(names=(char **)realloc(names, (i+1)*sizeof(char *))))
      errmess("realloc(names)");
    if (!(names[i]=(char *)calloc(strlen(line)+1, sizeof(char))))
      errmess("calloc(names[i])");

    sscanf(line, "%s", names[i]);
  }

  fclose(inf);

  *imnames=names;

  return(i);
}
/*--------------------------------------------------------*/
float *trim_image(int *pnx, int *pny, float **data, int verbose)
{
        int   nx0, ny0,
              nx, ny,
              i, j,
              x1, y1;
        float *nldata,
              *ldata;

  nx0=*pnx;
  ny0=*pny;
  ldata=*data;

  nx=(int)rint(pow(2.0, floor(log(nx0)/log(2.0))));
  ny=(int)rint(pow(2.0, floor(log(ny0)/log(2.0))));

  x1=(int)((nx0-nx)/2.0);
  y1=(int)((ny0-ny)/2.0);

  if (verbose > 1)
  {
    printf("nx: %d -> %d\n", nx0, nx);
    printf("ny: %d -> %d\n", ny0, ny);
    printf("x1= %d  y1= %d\n", x1, y1);
  }

  if (!(nldata=(float *)calloc(nx*ny, sizeof(float))))
    errmess("trim_image: calloc(nldata)");

  for (i=0; i<ny; i++)
    for (j=0; j<nx; j++)
      nldata[i*nx+j]=ldata[(i+y1)*nx0+j+x1];

  free(ldata);

  *pnx=nx;
  *pny=ny;

  return(nldata);
}
/*--------------------------------------------------------*/
int main(int argc, char *argv[])
{
        char    *imname1,     /* template image's file name   */
                *imname2,     /* list of images file name     */
                *instrname,   /* instrument parameters file   */
                *parfname,    /* parameters file name         */
                *outfname,    /* shifts file name             */
                **imnames,    /* input image names            */
                **header,     /* FITS header buffer           */
                outcfname[NAME_LEN];  /* cross-corelation map file name */
        int     nx0, ny0,     /* image dimensions             */
                i, j,         /* loop numerators              */
                n,
                tnn[2], cnn[2],
                nim,          /* number of images             */
                headlen,      /* number of FITS header lines  */
                csize,        /* size of cdata matrices       */
                tsize,        /* half of the tdatas' size     */
                corsize;      /* correlation tables size      */
        float   *cutim1,      /* template sector's buffer     */
                *cutim2,      /* image sector's buffer        */
                *tdata1,      /* template's FFT buffer        */
                *tdata2,      /* images' FFT buffer           */
                *corim8,      /* correlation matrix ?         */
                *corim4,      /* correlation matrix ?         */
                *cdata1,      /* template sector's FFT buffer */
                *cdata2,      /* image sector's FFT buffer    */
                x_shift,      /* x-shift between images       */
                y_shift,      /* y-shift between images       */
                x_max,        /* x-coord. of max. correlation */
                y_max,        /* y-coord. of max. correlation */
                rat,          /* see crf() */
                *data;        /* image data buffer            */
        FILE    *outf,        /* output shifts file descriptor*/
                *outcf;       /* correlation file descriptor  */
        PARAMS  par, cpar;

/* IO stuff */

  if ((argc < 6) || (argc > 7)) usage();

  parfname = argv[1];
  instrname = argv[2];
  outfname = argv[3];
  imname1  = argv[4];

  get_params(parfname, instrname, &par);
  if (par.verbose > 1)
  {
    printf("instrname: %s\n", instrname);
    printf("outfname: %s\n", outfname);
    printf("imname1:  %s\n", imname1);
  }

  if (!strcmp(argv[5], "-f"))
  {
    imname2=argv[6];
    if (par.verbose) printf("%s:\n", imname2);
    nim=readlist(imname2, &imnames);
  }
  else
  {
    if (!(imnames=(char **)calloc(1, sizeof(char *))))
      errmess("calloc(imnames)");
    if (!(imnames[0]=(char *)calloc(strlen(argv[5])+1, sizeof(char))))
      errmess("calloc(imnames[0])");
    strcpy(imnames[0], argv[5]);
    nim = 1;
  }

  if (par.verbose)
  {
    printf("%d images\n", nim);
    if (par.verbose > 1)
      for (i=0; i<nim; i++) printf("%s\n", imnames[i]);
    printf("Template: %s\n", imname1);
  }

  data=read_FITS_2D1file(imname1, 's', &headlen, &header, &par.nx0, &par.ny0);
  for (j=0; j<headlen; j++) free(header[j]);
  free(header);

  n=pow(2.0, floor(log(par.nx0)/log(2.0)));
  if (n != par.nx0)
  {
    data=trim_image(&par.nx0, &par.ny0, &data, par.verbose);
    if (par.verbose)
      printf("Template image trimmed to: %d x %d\n", par.nx0, par.ny0);
  }
  else
  {
    n=pow(2.0, floor(log(par.ny0)/log(2.0)));
    if (n != par.ny0)
    {
      data=trim_image(&par.nx0, &par.ny0, &data, par.verbose);
      if (par.verbose)
        printf("Template image trimmed to: %d x %d\n", par.nx0, par.ny0);
    }
  }

/* set parameters */
  cpar.sat_level = par.sat_level;
  cpar.min_level = par.min_level;
  cpar.bkg_frac  = par.bkg_frac;

  tnn[1] = par.nx = par.nx0/par.x_nbin;
  tnn[0] = par.ny = par.ny0/par.y_nbin;

  if (par.nx_cut > par.nx0) par.nx_cut=par.nx0;
  if (par.ny_cut > par.ny0) par.ny_cut=par.ny0;

  cnn[1] = cpar.nx = cpar.nx0 = par.nx_cut;
  cnn[0] = cpar.ny = cpar.nx0 = par.ny_cut;

  if (par.verbose > 1)
  {
    printf("tnn[1]= %d   tnn[0]= %d\n", tnn[1], tnn[0]);
    printf("cnn[1]= %d   cnn[0]= %d\n", cnn[1], cnn[0]);
  }

  cpar.x_nbin = 1;
  cpar.y_nbin = 1;

  par.dx = par.nx0/2 - par.nx_cut/2;
  par.dy = par.ny0/2 - par.ny_cut/2;

/* get memory */
  tsize=par.nx*par.ny;
  csize=cpar.nx*cpar.ny;
  corsize= (tsize > csize) ? tsize : csize;

  if (!(tdata1=(float *)malloc(2*tsize*sizeof(float))))
    errmess("malloc(tdata1)");
  if (!(tdata2=(float *)malloc(2*tsize*sizeof(float))))
    errmess("malloc(tdata2)");

  if (!(corim8=(float *)malloc(2*corsize*sizeof(float))))
    errmess("malloc(corim8)");
  if (!(corim4=(float *)malloc(  corsize*sizeof(float))))
    errmess("malloc(corim4)");

  if (!(cutim1=(float *)malloc(csize*sizeof(float))))
    errmess("malloc(cutim1)");
  if (!(cutim2=(float *)malloc(csize*sizeof(float))))
    errmess("malloc(cutim2)");

  if (!(cdata1=(float *)malloc(2*csize*sizeof(float))))
    errmess("malloc(cdata1)");
  if (!(cdata2=(float *)malloc(2*csize*sizeof(float))))
    errmess("malloc(cdata2)");

/* prepare stuff for template (once for all) */
  pack_data(data, tdata1, par);

  fourn(tdata1, tnn, 2, 1);

  get_sector(data, cutim1, par);

  pack_data(cutim1, cdata1, cpar);
  free(cutim1);

  fourn(cdata1, cnn, 2, 1);

  free(data);

  if (!(outf=fopen(outfname, "w")))  errmess(outfname);

/* loop over test images */
  for(i=0; i<nim; i++)
  {
/* get test image */
    if (par.verbose) printf("Image:    %s\n", imnames[i]);
    data=read_FITS_2D1file(imnames[i], 's', &headlen, &header, &nx0, &ny0);

    n=pow(2.0, floor(log(nx0)/log(2.0)));
    if (n != nx0)
    {
      data=trim_image(&nx0, &ny0, &data, par.verbose);
      if (par.verbose)
        printf("Image trimmed to: %d x %d\n", nx0, ny0);
    }
    else
    {
      n=pow(2.0, floor(log(ny0)/log(2.0)));
      if (n != ny0)
      {
        data=trim_image(&nx0, &ny0, &data, par.verbose);
        if (par.verbose)
          printf("Image trimmed to: %d x %d\n", nx0, ny0);
      }
    }

    for (j=0; j<headlen; j++) free(header[j]);
    free(header);

    if ((par.nx0 != nx0) || (par.ny0 != ny0))
    {
      printf("ERROR! Image size mismatch\n");
      if (par.verbose)
      {
        printf("  %s:  %d x %d\n", imname1, par.nx0, par.ny0);
        printf("  %s:  %d x %d\n", imnames[i], nx0, ny0);
      }
      exit(2);
    }

    pack_data(data, tdata2, par);

    rat = crf(tdata1, tdata2, corim8, tnn, 0);

    get_max(corim8, corim4, par.nx, par.ny, &x_max, &y_max);

    x_shift = (x_max - par.nx/2 )*par.x_nbin;
    y_shift = (y_max - par.ny/2 )*par.y_nbin;

    if (par.verbose) printf("first guess: dx= %f dy= %f\n", x_shift, y_shift);

/* if no binning first guess is the final answer */

    if ((par.x_nbin != 1) || (par.y_nbin != 1))
    {
      par.dx = par.nx0/2 - cpar.nx/2 - (int)x_shift;
      par.dy = par.ny0/2 - cpar.ny/2 - (int)y_shift;

      get_sector(data, cutim2, par);

      pack_data(cutim2, cdata2, cpar);

      rat = crf(cdata1, cdata2, corim8, cnn, 0);

      get_max(corim8, corim4, cpar.nx, cpar.ny, &x_max, &y_max);

      x_shift = x_max - cpar.nx/2 + (int)x_shift;
      y_shift = y_max - cpar.ny/2 + (int)y_shift;

      if (par.verbose)
        printf("final guess: dx= %f dy= %f\n", x_shift, y_shift);

      x_shift = floor(x_max - cpar.nx/2 + (int)x_shift + 0.5);
      y_shift = floor(y_max - cpar.ny/2 + (int)y_shift + 0.5);

      fprintf(outf, "%s  %4d  %4d\n", imnames[i], (int)x_shift, (int)y_shift);

      if (par.verbose > 2)
      {
        sprintf(outcfname, "%s.bmp", imnames[i]);
        if (!(outcf=fopen(outcfname, "w"))) errmess(outcfname);
        fwrite(corim4, csize*sizeof(float), 1, outcf);
        fclose(outcf);
      }
    }
    else
    {
      if (par.verbose)
        printf("final guess: dx= %f dy= %f\n", x_shift, y_shift);
      fprintf(outf, "%s  %4d  %4d\n", imnames[i], (int)x_shift, (int)y_shift);
    }

    free(data);
  }

  fclose(outf);

  for (i=0; i<nim; i++) free(imnames[i]);
  free(imnames);

  free(cutim2);
  free(tdata1);
  free(tdata2);
  free(corim8);
  free(corim4);
  free(cdata1);
  free(cdata2);

  return(0);
}
/*** END ***/
