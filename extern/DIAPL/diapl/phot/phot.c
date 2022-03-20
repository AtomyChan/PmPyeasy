/*========================================================*/
/*                                                        */
/*  phot.c              version 1.8.1   2006.01.02        */
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
/* Program does psf photometry on difference images.      */
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
#include <math.h>
#include <stdlib.h>

#include "defs.h"
#include "funcs.h"

#include "errmess.h"
#include "pfitsin1.h"

#define INIT_NIM 512

/*--------------------------------------------------------*/
int readfnames(char *iname, char ***diffiles, char ***imfiles, char ***kerfiles)
{
        char  **dif,
              **imf,
              **kef,
              line[501],
              s1[201], s2[201], s3[201];
        int   i;
        FILE  *inf;

  if (!(dif=(char **)malloc(sizeof(char *))))
    errmess("malloc(dif)");
  if (!(imf=(char **)malloc(sizeof(char *))))
    errmess("malloc(imf)");
  if (!(kef=(char **)malloc(sizeof(char *))))
    errmess("malloc(kef)");

  if (!(inf = fopen(iname, "r"))) errmess(iname);

  for (i=0; !feof(inf); i++)
  {
    fgets(line, 500, inf);
    if (feof(inf)) break;

    sscanf(line, "%s %s %s", s1, s2, s3);

    if (!(dif=(char **)realloc(dif, (i+1)*sizeof(char *))))
      errmess("realloc(dif)");
    if (!(dif[i]=(char *)calloc(strlen(s1)+1, sizeof(char))))
      errmess("calloc(dif[i])");

    if (!(imf=(char **)realloc(imf, (i+1)*sizeof(char *))))
      errmess("realloc(imf)");
    if (!(imf[i]=(char *)calloc(strlen(s2)+1, sizeof(char))))
      errmess("calloc(imf[i])");

    if (!(kef=(char **)realloc(kef, (i+1)*sizeof(char *))))
      errmess("realloc(kef)");
    if (!(kef[i]=(char *)calloc(strlen(s3)+1, sizeof(char))))
      errmess("calloc(kef[i])");

    strcpy(dif[i], s1);
    strcpy(imf[i], s2);
    strcpy(kef[i], s3);
  }

  fclose(inf);

  *diffiles=dif;
  *imfiles=imf;
  *kerfiles=kef;

  return(i);
}
/*--------------------------------------------------------*/
int readcoo(char *iname, int **vnum, double **x, double **y)
{
        char    line[201];
        int     i,
                *lvnum;
        double  *lx, *ly;
        FILE    *inf;

  if (!(lvnum=(int *)malloc(sizeof(int)))) errmess("malloc(lvnum)");
  if (!(lx=(double *)malloc(sizeof(double)))) errmess("malloc(lx)");
  if (!(ly=(double *)malloc(sizeof(double)))) errmess("malloc(ly)");

  if (!(inf = fopen(iname, "r"))) errmess(iname);

  for (i=0; !feof(inf); i++)
  {
    fgets(line, 200, inf);
    if (feof(inf)) break;

    if (!(lx=(double *)realloc(lx, (i+1)*sizeof(double))))
      errmess("realloc(lx)");
    if (!(ly=(double *)realloc(ly, (i+1)*sizeof(double))))
      errmess("realloc(ly)");
    if (!(lvnum=(int *)realloc(lvnum, (i+1)*sizeof(int))))
      errmess("realloc(vnum)");

    sscanf(line, "%d %lf %lf", &lvnum[i], &lx[i], &ly[i]);
  }

  fclose(inf);

  *vnum=lvnum;
  *x=lx;
  *y=ly;

  return(i);
}
/*--------------------------------------------------------*/
void writedba(char *oname, char **diffiles, int nim, int npsf, int *vnum,
              STAR **obj)
{
        int   i, j;
        FILE  *outf;
        STAR  *objp;

  if (!(outf=fopen(oname, "w"))) errmess(oname);

  fprintf(outf, "image     var_num  num_bad_pix  profile_flux  flux_err");
  fprintf(outf, "  ap_flux  flux_err  bkg  psf_chi^2(s)          fwhm\n");

  for (j=0; j<nim; j++)
  {
    for (i=0; i<npsf; i++)
    {
      objp = &obj[i][j];

      fprintf(outf, "%s  %5d  %3d  ", diffiles[j], vnum[i], objp->nbad);
      fprintf(outf, "%9g  %8g  ", objp->p_flux, objp->p_err);
      fprintf(outf, "%9g  %8g  ", objp->a_flux, objp->a_err);
      fprintf(outf, "%8g  %8g  %g  ", objp->bg, objp->chi2_n, objp->corr);
      fprintf(outf, "%9g  %8g\n", objp->ker_chi2, objp->fwhm);
    }
  }

  fclose(outf);

  return;
}
/*--------------------------------------------------------*/
void writedbb(char *oname, char **diffiles, int nim, int npsf, int *vnum,
              STAR **obj)
{
        char  record[10*sizeof(double)];
        int   i, j,
              rs,
              ofs;
        FILE  *outf;
        STAR  *objp;

  rs=sizeof(double);

  if (!(outf=fopen(oname, "w"))) errmess(oname);

  for(j=0; j<nim; j++)
  {
    for (i=0; i<npsf; i++)
    {
      objp = &obj[i][j];

      ofs=0;
      memcpy(&record[ofs], &objp->p_flux, rs);    ofs+=rs;
      memcpy(&record[ofs], &objp->p_err,  rs);    ofs+=rs;
      memcpy(&record[ofs], &objp->a_flux, rs);    ofs+=rs;
      memcpy(&record[ofs], &objp->a_err,  rs);    ofs+=rs;
      memcpy(&record[ofs], &objp->bg,     rs);    ofs+=rs;
      memcpy(&record[ofs], &objp->chi2_n, rs);    ofs+=rs;
      memcpy(&record[ofs], &objp->corr,   rs);    ofs+=rs;
      memcpy(&record[ofs], &objp->ker_chi2, rs);  ofs+=rs;
      memcpy(&record[ofs], &objp->fwhm,   rs);    ofs+=rs;
      memcpy(&record[ofs], &objp->nbad,   sizeof(int));

      fwrite(record, 1, sizeof(record), outf);
    }
  }

  fclose(outf);

  return;
}
/*--------------------------------------------------------*/
int main(int argc, char *argv[])
{
        char        **newheader,
                    *imffname, *parfname, *psffname, *reffname,
                    *instrname,
                    *fieldname,
                    *outfname,
                    *coofname,
                    **diffiles,
                    **imfiles,
                    **kerfiles;
        int         nx0, ny0, k, *sector, nim, iim, npsf, ipsf, hsize,
                    isec_x, isec_y, *cx, *cy, psfn, kern, irad,
//                  ofs,
                    *vnum,
                    i;
        double       *im, *difim, *refim;
        double      **wxy, *x, *y, *xs, *ys, **psfs, *psfim, *kerim, ratio;
        STAR        **obj, *objp;
        PSF_STRUCT  psf;
        KER_STRUCT  ker;
        PAR_STRUCT  par;

/*** IO stuff ***/

  if (argc != 7)
  {
    printf("\n\tUSAGE: phot  parameter_file instrument_file");
    printf(" ref_image psf_fit_file image_data_list field_name\n");
    exit(1);
  }

  parfname= argv[1];
  instrname=argv[2];
  reffname= argv[3];
  psffname= argv[4];
  imffname= argv[5];
  fieldname=argv[6];

  get_params(parfname, instrname, &par);

  if (par.verbose)
    printf("\n\n*** Profile photometry with variable PSF and kernel ***\n\n");

  if (!(outfname=(char *)calloc(strlen(fieldname)+5, sizeof(char))))
    errmess("calloc(outfname)");
  strcpy(outfname, fieldname);
  strcat(outfname, ".db");

  if (!(coofname=(char *)calloc(strlen(fieldname)+5, sizeof(char))))
    errmess("calloc(coofname)");
  strcpy(coofname, fieldname);
  strcat(coofname, ".coo");

/***  read filenames for difference images, images and kernel fits  ***/
  nim=readfnames(imffname, &diffiles, &imfiles, &kerfiles);
  if (par.verbose)  printf("%d images to process\n", nim);

/***  read coordinates of variables  ***/
  npsf=readcoo(coofname, &vnum, &x, &y);
  if (par.verbose)  printf("%d variables to measure\n\n", npsf);

/***  read in psf fit and get a sample kernel from the first image  ***/
  read_psf(psffname, &psf);
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

  if (!(sector=(int *)malloc(npsf*sizeof(int)))) errmess("malloc(sector)");
  if (!(cx= (int *)malloc(npsf*sizeof(int))))          errmess("malloc(cx)");
  if (!(cy= (int *)malloc(npsf*sizeof(int))))          errmess("malloc(cy)");
  if (!(xs= (double *)malloc(npsf*sizeof(double))))    errmess("malloc(xs)");
  if (!(ys= (double *)malloc(npsf*sizeof(double))))    errmess("malloc(ys)");
  if (!(wxy=(double **)malloc(npsf*sizeof(double *)))) errmess("malloc(wxy)");

  if (!(psfs=(double **)malloc(npsf*sizeof(double *))))
    errmess("malloc(psfs)");
  if (!(kerim=(double *)malloc(kern*kern*sizeof(double))))
    errmess("malloc(kerim)");
  if (!(psfim=(double *)malloc(psfn*psfn*sizeof(double))))
    errmess("malloc(psfim)");

  if (!(obj=(STAR **)malloc(npsf*sizeof(STAR *)))) errmess("malloc(obj)");

/***********************************************************************/
/** get things that can be done once for all: spatial coeffs and psfs **/
/***********************************************************************/

  for (ipsf=0; ipsf<npsf; ipsf++)
  {
    if (!(obj[ipsf]=(STAR *)malloc(nim*sizeof(STAR))))
      errmess("malloc(obj[ipsf])");
    if (!(wxy[ipsf]=(double *)malloc(ker.nwxy*sizeof(double))))
      errmess("malloc(wxy[ipsf])");
    if (!(psfs[ipsf]=(double *)malloc(psfn*psfn*sizeof(double))))
      errmess("malloc(psfs[ipsf])");

/* offsets for image sectors from kernel fit */

    cx[ipsf]=(int)floor(x[ipsf]+0.5);
    cy[ipsf]=(int)floor(y[ipsf]+0.5);

    isec_x=(cx[ipsf] - ker.hw)/(ker.nx - 2*ker.hw);
    isec_y=(cy[ipsf] - ker.hw)/(ker.ny - 2*ker.hw);

    xs[ipsf]=x[ipsf] - isec_x*(ker.nx - 2*ker.hw);
    ys[ipsf]=y[ipsf] - isec_y*(ker.ny - 2*ker.hw);

    sector[ipsf]=isec_y + ker.nsec_y*isec_x;

    spatial_coeffs(&ker, xs[ipsf], ys[ipsf], wxy[ipsf]);

    init_psf(&psf, x[ipsf], y[ipsf]);
    make_psf(&psf, x[ipsf], y[ipsf], cx[ipsf], cy[ipsf], psfs[ipsf]);
  }

  refim=read_FITS_2D1file(reffname, 's', &hsize, &newheader, &nx0, &ny0);
  for (i=0; i<hsize; i++) free(newheader[i]);
  free(newheader);

  make_vectors(&ker);

  par.nx0=nx0;
  par.ny0=ny0;
  par.psfhw=psf.hw;
  par.psfn=psfn;

  irad=(int)par.anrad2 + 2;

/*******************************/
/***  main loop over images  ***/
/*******************************/

  for (iim=0; iim<nim; iim++)
  {
    if (par.verbose > 2) printf("%d: %s\n", iim, diffiles[iim]);

    difim=read_FITS_2D1file(diffiles[iim], 's', &hsize, &newheader, &nx0, &ny0);

    for (i=0; i<hsize; i++) free(newheader[i]);
    free(newheader);

    if ((nx0 != par.nx0) || (ny0 != par.ny0))
    {
      printf("ERROR! phot: wrong size of the image %s\n", diffiles[iim]);
      exit(2);
    }

    im=read_FITS_2D1file(imfiles[iim], 's', &hsize, &newheader, &nx0, &ny0);
    for (i=0; i<hsize; i++) free(newheader[i]);
    free(newheader);

    if ((nx0 != par.nx0) || (ny0 != par.ny0))
    {
      printf("ERROR! phot: wrong size of the image %s\n", imfiles[iim]);
      exit(3);
    }

/* read kernel into tables allocated before */
    read_kernel(kerfiles[iim], &ker, 0, par.verbose);

/*** loop over stars ***/
    for (ipsf=0; ipsf<npsf; ipsf++)
    {
      objp = &obj[ipsf][iim];

      if ((cx[ipsf] < irad) || (cy[ipsf] < irad) ||
          (cx[ipsf] >= nx0-irad) || (cy[ipsf] >= ny0-irad))
      {
        if (par.verbose)
          printf("%s warning: object %4d too close to edge: ignored!\n",
                  diffiles[iim], ipsf);

        objp->a_flux   = par.bad_value;
        objp->a_err    = par.bad_value;
        objp->p_flux   = par.bad_value;
        objp->p_err    = par.bad_value;
        objp->chi2_n   = par.bad_value;
        objp->ker_chi2 = par.bad_value;
        objp->corr     = par.bad_value;
        objp->nbad     = -1;

        continue;
      }

/*** prepare local psf ***/
      make_kernel(&ker, wxy[ipsf], kerim, sector[ipsf]);

      im_convolve(psfs[ipsf], psfim, psfn, psfn, kerim, ker.hw);

      objp->ker_chi2=(double)ker.chi2_n[sector[ipsf]];

      objp->fwhm=(double)get_fwhm(psfim, x[ipsf], y[ipsf], cx[ipsf], cy[ipsf],
                                 &par, &ratio);

      if (par.bkg_mode) objp->bg = bkg(difim, cx[ipsf], cy[ipsf], &par);
      else              objp->bg = 0.0;

/*** actual profile and aperture photometry ***/
      get_phot(difim, im, refim, x[ipsf], y[ipsf], cx[ipsf], cy[ipsf], psfim,
               objp, &par);

      if (par.verbose > 1)
        printf("%s  star: %5d  flux= %9g +- %8g   nbad: %d\n",
                diffiles[iim], ipsf, objp->p_flux, objp->p_err, objp->nbad);
    }

    free(difim);
    free(im);
  }

/*** write photometry to the output file  ***/
  if (par.verbose) printf("\nWriting photometry to:  %s\n", outfname);

  if (par.dbf == 'A')
    writedba(outfname, diffiles, nim, npsf, vnum, obj);
  else
    writedbb(outfname, diffiles, nim, npsf, vnum, obj);

  if (par.verbose)  printf("\nPhotometry done!\n");

  return(0);
}
/*** END ***/
