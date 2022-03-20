/*========================================================*/
/*                                                        */
/*  getvar.c            version 1.3.1   2005.04.20        */
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
/*  Finds variables in a series of difference images.     */
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

#include "defs.h"
#include "funcs.h"

#include "errmess.h"
#include "pfitsio1.h"

#define RECORD_LEN 52

/*--------------------------------------------------------*/
void usage()
{
  printf("\n\tUSAGE: getvar  parameter_file instrument_file ref_image ");
  printf("psf_fit_file image_data_list field_name x0_off y0_off\n\n");
  exit(1);
}
/*--------------------------------------------------------*/
/*  Read filenames for:                                   */
/*       difference images, images and kernel fit params  */
/*  Return number of images                               */
/*--------------------------------------------------------*/
int read_inp_list(char *listname,
                  char ***difnames, char ***imnames, char ***kernames)
{
        char  **ldifnames,
              **limnames,
              **lkernames,
              line[501],
              s1[201], s2[201], s3[201];
        int   i;
        FILE  *inf;

  if (!(ldifnames=(char **)calloc(1, sizeof(char *))))
    errmess("calloc(ldifnames)");
  if (!(limnames=(char **)calloc(1, sizeof(char *))))
    errmess("calloc(limnames)");
  if (!(lkernames=(char **)calloc(1, sizeof(char *))))
    errmess("calloc(lkernames)");

  if (!(inf=fopen(listname, "r"))) errmess(listname);

  for (i=0; !feof(inf); i++)
  {
    fgets(line, 200, inf);
    if (feof(inf)) break;

    sscanf(line, "%s %s %s", s1, s2, s3);

    if (!(ldifnames=(char **)realloc(ldifnames, (i+1)*sizeof(char *))))
      errmess("realloc(ldifnames)");
    if (!(limnames=(char **)realloc(limnames, (i+1)*sizeof(char *))))
      errmess("realloc(limnames)");
    if (!(lkernames=(char **)realloc(lkernames, (i+1)*sizeof(char *))))
      errmess("realloc(lkernames)");

    if (!(ldifnames[i]=(char *)calloc(strlen(s1)+1, sizeof(char))))
      errmess("calloc(ldifnames[i])");
    if (!(limnames[i]=(char *)calloc(strlen(s2)+1, sizeof(char))))
      errmess("calloc(limnames[i])");
    if (!(lkernames[i]=(char *)calloc(strlen(s3)+1, sizeof(char))))
      errmess("calloc(lkernames[i])");

    strcpy(ldifnames[i], s1);
    strcpy(limnames[i], s2);
    strcpy(lkernames[i], s3);
  }

  fclose(inf);

  *difnames=ldifnames;
  *imnames=limnames;
  *kernames=lkernames;

  return(i);
}
/*--------------------------------------------------------*/
int make_header(char ***header, int nx, int ny)
{
        char **lheader, card[CARD_SIZE+1];
        int  i;

  if (!(lheader=(char **)calloc(6, sizeof(char *)))) errmess("calloc(lheader)");
  for (i=0; i<6; i++)
    if (!(lheader[i]=(char *)calloc(CARD_SIZE, sizeof(char))))
      errmess("calloc(lheader[i])");

  strcpy(card, "SIMPLE  =                    T  /  ");
  strcat(card, "FITS STANDARD                                ");
  memcpy(lheader[0], card, CARD_SIZE);

  strcpy(card, "BITPIX  =                  -32  /  ");
  strcat(card, "FITS BITS/PIXEL                              ");
  memcpy(lheader[1], card, CARD_SIZE);

  strcpy(card, "NAXIS   =                    2  /  ");
  strcat(card, "NUMBER OF AXES                               ");
  memcpy(lheader[2], card, CARD_SIZE);

  sprintf(card, "NAXIS1  =               %6d  /  ", nx);
  strcat(card, "                                             ");
  memcpy(lheader[3], card, CARD_SIZE);

  sprintf(card, "NAXIS2  =               %6d  /  ", ny);
  strcat(card, "                                             ");
  memcpy(lheader[4], card, CARD_SIZE);

  strcpy(card, "END                                ");
  strcat(card, "                                             ");
  memcpy(lheader[5], card, CARD_SIZE);

  *header=lheader;

  return(6);
}
/*--------------------------------------------------------*/
void writevar2fits(char *fname, int nx, int ny, float *data)
{
        char  **header;
        int   hsize, i;

  hsize=make_header(&header, nx, ny);

  write_FITS_2D1file(fname, hsize, header, nx*ny, sizeof(float), data);

  for (i=0; i<hsize; i++) free(header[i]);
  free(header);
  
  return;
}
/*--------------------------------------------------------*/
int main(int argc, char *argv[])
{
        char        **header,
                    *imffname, *parfname, *psffname, *reffname,
                    *instrname,
                    *fieldname,
                    *catname,
                    *coofname,
                    record[RECORD_LEN],
                    **diffiles,
                    **imfiles,
                    **kerfiles;
        int         nx0, ny0, k, nim, iim, hsize,
                    cx, cy, psfn, kern, irad, ofs, nobj1, nobj1_max, nobj2,
                    nobj2_max, i, j, flag, *index,
                    nvar;
        float       x0_off, y0_off, x_tmp, y_tmp, fwhm_limit, *fwhm,
                    **im, **difim, *refim, *varim1, *varim2, *corrim,
                    *tmpim;
        double      *wxy, x, y, *psfs, *psfim, *kerim, ratio;
        FILE        *outfcat, *outfcoo;
        STAR        *obj1, *obj2, *objp;
        PSF_STRUCT  psf;
        KER_STRUCT  ker;
        PAR_STRUCT  par;

/*** IO stuff ***/

  if (argc != 9) usage();

  parfname = argv[1];
  instrname= argv[2];
  reffname = argv[3];
  psffname = argv[4];
  imffname = argv[5];
  fieldname = argv[6]; 
  x0_off = atof(argv[7]);
  y0_off = atof(argv[8]);

  get_params(parfname, instrname, &par);

  if (!(catname=(char *)calloc(strlen(fieldname)+5, sizeof(char))))
    errmess("calloc(catname)");
  strcpy(catname, fieldname);   strcat(catname, ".cat");

  if (!(coofname=(char *)calloc(strlen(fieldname)+5, sizeof(char))))
    errmess("calloc(coofname)");
  strcpy(coofname, fieldname);  strcat(coofname, ".coo");

  if (par.verbose > 2)
  {
    printf("parfname = %s\n", parfname);
    printf("instrname= %s\n", instrname);
    printf("reffname = %s\n", reffname);
    printf("psffname = %s\n", psffname);
    printf("imffname = %s\n", imffname);
    printf("fieldname= %s\n", fieldname);
    printf("x0_off   = %g\n", x0_off);
    printf("y0_off   = %g\n", y0_off);
    printf("--------\n");
    printf("catname  = %s\n", catname);
    printf("coofname = %s\n", coofname);
    printf("--------\n");
  }

  nim=read_inp_list(imffname, &diffiles, &imfiles, &kerfiles);
  if (par.verbose)  printf("%s: %d images\n", imffname, nim);

/**********************************************************************/
/***  read in psf fit and get a sample kernel from the first image  ***/
/**********************************************************************/

  read_psf(psffname, &psf, par.verbose);
  read_kernel(kerfiles[0], &ker, 1, par.verbose);

  psf.normrad = par.normrad;
  psf.hw += ker.hw;

  psfn = 2*psf.hw + 1;
  kern = 2*ker.hw + 1;

/*** get memory ***/
  if (!(ker.vecs = (double **)malloc(ker.nvecs*sizeof(double *))))
    errmess("malloc(ker.vecs)");
  for (k=0; k<ker.nvecs; k++)
    if (!(ker.vecs[k] = (double *)malloc(kern*kern*sizeof(double))))
      errmess("malloc(ker.vecs[k])");

  if (!(kerim=(double *)malloc(kern*kern*sizeof(double))))
    errmess("malloc(kerim)");
  if (!(psfim=(double *)malloc(psfn*psfn*sizeof(double))))
    errmess("malloc(psfim)");

  if (!(fwhm  = (float *)malloc(nim*sizeof(float)))) errmess("malloc(fwhm)");
  if (!(index = (int   *)malloc(nim*sizeof(int)))) errmess("malloc(index)");

/***********************************************************************/
/** get things that can be done once for all: spatial coeffs and psfs **/
/***********************************************************************/

  refim=read_FITS_2D1file(reffname, 's', &hsize, &header, &nx0, &ny0);
  for (i=0; i<hsize; i++) free(header[i]);
  free(header);

  par.nx0 = nx0;
  par.ny0 = ny0;
  par.psfhw = psf.hw;
  par.psfn  = psfn;

  irad  = (int)par.anrad2 + 2;

/*** get even more memory ***/
  if (!(psfs=(double *)malloc(psfn*psfn*sizeof(double))))
    errmess("malloc(psfs)");
  if (!(wxy=(double *)malloc(ker.nwxy*sizeof(double)))) errmess("malloc(wxy)");

  if (!(im   =(float **)malloc(nim*sizeof(float *)))) errmess("malloc(im)");
  if (!(difim=(float **)malloc(nim*sizeof(float *)))) errmess("malloc(difim)");

  if (!(varim1=(float *)malloc(nx0*ny0*sizeof(float))))
    errmess("malloc(varim1)");
  if (!(varim2=(float *)malloc(nx0*ny0*sizeof(float))))
    errmess("malloc(varim2)");
  if (!(corrim=(float *)malloc(nx0*ny0*sizeof(float))))
    errmess("malloc(corrim)");

  for (i=0; i<nx0*ny0; i++) varim1[i] = varim2[i] = 0.0;

/*** make reference psf and spatial coeffs for kernel ***/
  init_psf(&psf, (double)(nx0/2), (double)(ny0/2));
  make_psf(&psf, (double)(nx0/2), (double)(ny0/2), nx0/2, ny0/2, psfs);

  make_vectors(&ker);

  spatial_coeffs(&ker, (double)(ker.nx/2), (double)(ker.ny/2), wxy);

/*******************************/
/***  main loop over images  ***/
/*******************************/

  for (iim=0; iim<nim; iim++)
  {
    read_kernel(kerfiles[iim], &ker, 0, par.verbose);
    make_kernel(&ker, wxy, kerim, 0);
    im_convolve(psfs, psfim, psfn, psfn, kerim, ker.hw);
    fwhm[iim] = get_fwhm(psfim, 0.0, 0.0, 0, 0, &par, &ratio);
  }

/*** reject bad seeing frames ***/
  quick_sort(fwhm, index, nim);

  fwhm_limit = fwhm[index[(int)(nim*par.fwhm_frac)]];

  j = 0;
  for (iim=0; iim<nim; iim++)
  {
    if (fwhm[iim] <= fwhm_limit)
    {
      if (j != iim)
      {
        strcpy(diffiles[j], diffiles[iim]);
        strcpy( imfiles[j], imfiles[iim]);
        strcpy(kerfiles[j], kerfiles[iim]);
      }
      j++;
    }
  }

  nim = j;

  if (nim == 0)
  {
    printf("no images left ! (exit)\n");
    exit(2);
  }

/*** get all variability data at once ***/

  for (iim=0; iim<nim; iim++)
  {
    difim[iim]=read_FITS_2D1file(diffiles[iim], 's',
                                  &hsize, &header, &nx0, &ny0);
    for (i=0; i<hsize; i++) free(header[i]);
    free(header);

    if ((nx0 != par.nx0) || (ny0 != par.ny0))
    {
      printf("ERROR! getvar: image %s has wrong size\n", diffiles[iim]);
      exit(3);
    }

    im[iim]=read_FITS_2D1file(imfiles[iim], 's', &hsize, &header, &nx0, &ny0);
    for (i=0; i<hsize; i++) free(header[i]);
    free(header);

    if ((nx0 != par.nx0) || (ny0 != par.ny0))
    {
      printf("ERROR! getvar: image %s has wrong size\n", imfiles[iim]);
      exit(4);
    }
  }

/*** prepare variability likelihood image ***/
  get_repeaters(difim, im, varim1, varim2, &par, nim);

/*** save variability images ***/
  if (!(tmpim=(float *)malloc(nx0*ny0*sizeof(float)))) errmess("malloc(tmpim)");

  memcpy(tmpim, varim1, nx0*ny0*sizeof(float));
  writevar2fits("var1.fits", nx0, ny0, tmpim);

  memcpy(tmpim, varim2, nx0*ny0*sizeof(float));
  writevar2fits("var2.fits", nx0, ny0, tmpim);

  free(tmpim);

/*** find stellar looking things ***/
  covar_sig(varim1, corrim, nx0, ny0, 0.0, par.mohw, psfs, par.psfhw);

  obj1 = find_stars(corrim, nx0, ny0, &nobj1, &nobj1_max, &par);
  if (par.verbose > 1) printf("find_stars -> %d variables of type 1\n", nobj1);

  center_stars(obj1, nobj1, difim, im, &par, nim);

/***/
  covar_sig(varim2, corrim, nx0, ny0, 0.0, par.mohw, psfs, par.psfhw);

  obj2 = find_stars(corrim, nx0, ny0, &nobj2, &nobj2_max, &par);
  if (par.verbose > 1) printf("find_stars -> %d variables of type 2\n", nobj2);

  center_stars(obj2, nobj2, difim, im, &par, nim);

/*** some stars may be found in both: take only one type and flag ***/
  cross_id(obj1, nobj1, obj2, nobj2, par.id_rad);

  if (par.verbose)
    printf("\nFound %d and %d candidates for variables of both types\n\n",
     nobj1, nobj2);

/*********************************************************************/
/***   write results to a binary file and x,y to a temporary file  ***/
/*********************************************************************/

  if (!(outfcoo = fopen(coofname, "w"))) errmess(coofname);
  if (!(outfcat = fopen(catname, "w"))) errmess(catname);

  free(coofname);
  free(catname);

  fseek(outfcat, 0, SEEK_END);

  ofs = sizeof(float);

/*** get sinusoidal variables out ***/
  if (par.verbose) printf("var1:\n num     X       Y      flux      flag\n");
  nvar=0;

  for (i=0; i<nobj1; i++)
  {
    if (par.verbose > 1) printf("i= %d\n", i);

    objp = &(obj1[i]);

    x = (double)(objp->x);
    y = (double)(objp->y);
    if (par.verbose > 1) printf("x= %g   y= %g\n", x, y);
    if ((x < par.bad_margin) || (x >= par.nx0 - par.bad_margin)) continue;
    if ((y < par.bad_margin) || (y >= par.ny0 - par.bad_margin)) continue;

    cx = objp->cx;
    cy = objp->cy;
    if (par.verbose > 1) printf("cx= %d   cy= %d\n", cx, cy);

/*** append to a binary file ***/
    if (objp->nframes > 0)
    {
      fprintf(outfcoo, "%5d  %11.5f  %11.5f\n", nvar, x, y);
      nvar++;

/*** get template photometry ***/
      init_psf(&psf, x, y);

      make_psf(&psf, x, y, cx, cy, psfs);

      objp->bg = bkg(refim, cx, cy, &par);

      get_phot(refim, refim, refim, x, y, cx, cy, psfs, objp, &par);

      flag = neighbor(refim, (int)x, (int)y, &par, objp->bg);

      if (par.verbose)
        printf("%4d  %9.3f  %9.3f  %9.1f  %5d\n", i, x, y, objp->p_flux, flag);

      x_tmp = (float)(x + x0_off);
      y_tmp = (float)(y + y0_off);

      if (par.verbose > 2)
      {
        printf("x_tmp= %g   y_tmp= %g\n", x_tmp, y_tmp);
        printf("ofs= %d   13*ofs= %d\n", ofs, 13*ofs);
      }

      memcpy(&record[ 0*ofs], &x_tmp        , ofs);
      memcpy(&record[ 1*ofs], &y_tmp        , ofs);
      memcpy(&record[ 2*ofs], &objp->p_flux , ofs);
      memcpy(&record[ 3*ofs], &objp->p_err  , ofs);
      memcpy(&record[ 4*ofs], &objp->a_flux , ofs);
      memcpy(&record[ 5*ofs], &objp->a_err  , ofs);
      memcpy(&record[ 6*ofs], &objp->bg     , ofs);
      memcpy(&record[ 7*ofs], &objp->chi2_n , ofs);
      memcpy(&record[ 8*ofs], &objp->corr   , ofs);
      memcpy(&record[ 9*ofs], &objp->nbad   , sizeof(int));
      memcpy(&record[10*ofs], &objp->vtype  , sizeof(int));
      memcpy(&record[11*ofs], &objp->nframes, sizeof(int));
      memcpy(&record[12*ofs], &flag         , sizeof(int));

      if (par.verbose > 2) printf("record prepared - writting...\n");
      fwrite(record, 1, sizeof(record), outfcat);
    }
  }

/*** the same for transients ***/
  if (par.verbose) printf("var2:\n num     X       Y      flux      flag\n");

  for (i=0; i<nobj2; i++)
  {
    if (par.verbose > 1) printf("i= %d\n", i);

    objp = &obj2[i];

    x = (double)(objp->x);
    y = (double)(objp->y);
    if (par.verbose > 1) printf("x= %g   y= %g\n", x, y);
    if ((x < par.bad_margin) || (x >= par.nx0 - par.bad_margin)) continue;
    if ((y < par.bad_margin) || (y >= par.ny0 - par.bad_margin)) continue;

    cx = objp->cx;
    cy = objp->cy;
    if (par.verbose > 1) printf("cx= %d   cy= %d\n", cx, cy);

/*** append to a binary file ***/
    if (objp->nframes > 0)
    {
      fprintf(outfcoo, "%5d  %11.5f  %11.5f\n", nvar, x, y);
      nvar++;

/*** get template photometry ***/
      init_psf(&psf, x, y);

      make_psf(&psf, x, y, cx, cy, psfs);

      obj2[i].bg = bkg(refim, cx, cy, &par);

      get_phot(refim, refim, refim, x, y, cx, cy, psfs, objp, &par);

      flag = neighbor(refim, (int)x, (int)y, &par, objp->bg);

      if (par.verbose)
        printf("%4d  %9.3f  %9.3f  %9.1f  %5d\n",
                i+nobj1, x, y, objp->p_flux, flag);

      x_tmp = (float)(x + x0_off);
      y_tmp = (float)(y + y0_off);

      if (par.verbose > 2)
      {
        printf("x_tmp= %g   y_tmp= %g\n", x_tmp, y_tmp);
        printf("ofs= %d   13*ofs= %d\n", ofs, 13*ofs);
      }

      memcpy(&record[ 0*ofs], &x_tmp        , ofs);
      memcpy(&record[ 1*ofs], &y_tmp        , ofs);
      memcpy(&record[ 2*ofs], &objp->p_flux , ofs);
      memcpy(&record[ 3*ofs], &objp->p_err  , ofs);
      memcpy(&record[ 4*ofs], &objp->a_flux , ofs);
      memcpy(&record[ 5*ofs], &objp->a_err  , ofs);
      memcpy(&record[ 6*ofs], &objp->bg     , ofs);
      memcpy(&record[ 7*ofs], &objp->chi2_n , ofs);
      memcpy(&record[ 8*ofs], &objp->corr   , ofs);
      memcpy(&record[ 9*ofs], &objp->nbad   , sizeof(int));
      memcpy(&record[10*ofs], &objp->vtype  , sizeof(int));
      memcpy(&record[11*ofs], &objp->nframes, sizeof(int));
      memcpy(&record[12*ofs], &flag         , sizeof(int));

      if (par.verbose > 2) printf("record prepared - writting...\n");
      fwrite(record, 1, sizeof(record), outfcat);
    }
  }

  fclose(outfcoo);
  fclose(outfcat);

  if (par.verbose)  printf("\nVariability seach fihished!\n");

  return(0);
}
/*** END ***/
