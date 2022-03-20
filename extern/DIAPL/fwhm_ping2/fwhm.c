/*========================================================*/
/*                                                        */
/*  fwhm.c          version 2.8.5       2006.10.23        */
/*                                                        */
/*  Copyright (C) 2006 by Wojtek Pych, CAMK PAN           */
/*                                                        */
/*  Written for GNU project C and C++ Compiler            */
/*                                                        */
/*  Find FWHM for a given image.                          */
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

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>

#include "errmess.h"
#include "pfitsin.h"

/* input parameters structure */
typedef struct
{
  char  sky_mod[9],
        datatype[9];
  short vopt;
  int   aperture,
        iskyr,
        oskyr,
        margin;
  double minfwhm, maxfwhm,
        minpeak, maxpeak,
        minsky, maxsky;
} PARAMS;

void    usage();
PARAMS  readpar(char *parfname);
void    minmax(short, double **, int, int, double *, double *);
int     calc_mean(double *, int, double *, double *);
int     calc_mean2(double, double, double **, int, int, double *, double *);
void    make_hist(double **, int, int, double, double, int, double, int *);
double  calc_median(int *, int, int);
double  qlsq(size_t, double *);
double  calc_mode(PARAMS, size_t, int *);
double  stats(PARAMS, int, int, double **);
int     reject(short, double, double **, int, int *, int *);
int     find(PARAMS, double **, double, int, int, int **, int **);
void    centroid(PARAMS, double **, int, int, int, int, double, double *, double *);
char    plsq(short, double *, double *, double *, int, double *);
double  find_local_sky(PARAMS, double **, int, int, int, int);
char    fitgauss(PARAMS, double **, double, int, int, int, double, double,
                 double *, int);
void    sort(double *, int);
int     local_skies(PARAMS, double **, int, int, int, int **, int **, double **);
int     det_cent(PARAMS, double **, int, int, int, int *, int *, double *,
                 double **, double **);
int     fwhms(PARAMS, double **, int, int, int, double *, double *, double *,
              double *, double *, char *);

/*--------------------------------------------------------*/
void p_errmess(char *errmess)
{
  printf("ERROR! %s\n", errmess);
  exit(EXIT_FAILURE);
}
/*--------------------------------------------------------*/
void usage()
{
  printf("\n\t USAGE: fwhm parameter_file <frame> output_file\n\n");
  printf("\t output_file contain the star information from which the fwhm is derived\n");
  exit(EXIT_FAILURE);
}
/*--------------------------------------------------------*/
PARAMS readpar(char *parfname)
{
#define LLEN 200

        char    line[LLEN],
                key[LLEN],
                val[LLEN];
        int     i;
        FILE    *inf;
        PARAMS  lpar;

  if (!(inf = fopen(parfname, "r"))) errmess(parfname);

  for (i=0; !feof(inf); i++)
  {
    if (fgets(line, LLEN, inf) == NULL)
    {
      printf("ERROR! reading %s\n", parfname);
      exit(EXIT_FAILURE);
    }
    if (feof(inf)) break;

    sscanf(line, "%s = %s", key, val);

    if (!strcmp(key, "END")) break;
    else if (!strcasecmp(key, "SKY_MOD"))     sscanf(val, "%s", lpar.sky_mod);
    else if (!strcasecmp(key, "DATA_SIGN"))   sscanf(val, "%s", lpar.datatype);
    else if (!strcasecmp(key, "APERTURE"))    sscanf(val, "%d", &lpar.aperture);
    else if (!strcasecmp(key, "IN_SKY_RAD"))  sscanf(val, "%d", &lpar.iskyr);
    else if (!strcasecmp(key, "OUT_SKY_RAD")) sscanf(val, "%d", &lpar.oskyr);
    else if (!strcasecmp(key, "MARGIN"))      sscanf(val, "%d", &lpar.margin);
    else if (!strcasecmp(key, "MIN_FWHM"))    sscanf(val, "%lf", &lpar.minfwhm);
    else if (!strcasecmp(key, "MAX_FWHM"))    sscanf(val, "%lf", &lpar.maxfwhm);
    else if (!strcasecmp(key, "MIN_PEAK"))    sscanf(val, "%lf", &lpar.minpeak);
    else if (!strcasecmp(key, "MAX_PEAK"))    sscanf(val, "%lf", &lpar.maxpeak);
    else if (!strcasecmp(key, "MIN_SKY"))     sscanf(val, "%lf", &lpar.minsky);
    else if (!strcasecmp(key, "MAX_SKY"))     sscanf(val, "%lf", &lpar.maxsky);
    else if (!strcasecmp(key, "VERBOSE"))     sscanf(val, "%hd", &lpar.vopt);
    else printf("WARNING: Unknown parameter %s\n", key);
  }

  fclose(inf);

  if (lpar.vopt > 1)
  {
    printf("Input Parameters:\n");
    printf("-----------------\n");
    printf("Sky algorithm    = %s\n", lpar.sky_mod);
    printf("Data type        = %s\n", lpar.datatype);
    printf("Aperture size    = %d\n", lpar.aperture);
    printf("Inner Sky Radius = %d\n", lpar.iskyr);
    printf("Outer Sky Radius = %d\n", lpar.oskyr);
    printf("Margin           = %d\n", lpar.margin);
    printf("Minimum FWHM     = %lf\n", lpar.minfwhm);
    printf("Maximum FWHM     = %lf\n", lpar.maxfwhm);
    printf("Minimum peak     = %lf\n", lpar.minpeak);
    printf("Maximum peak     = %lf\n", lpar.maxpeak);
    printf("Minimum sky      = %lf\n", lpar.minsky);
    printf("Maximum sky      = %lf\n", lpar.maxsky);
    printf("Verbose          = %d\n", lpar.vopt);
    printf("-----------------\n");
  }

#undef LLLEN

  return(lpar);
}
/*--------------------------------------------------------*/
void minmax(short vopt, double **data, int xsize, int ysize,
            double *min, double *max)
{
        int i,                  /* rows loop numerator    */
            j,                  /* columns loop numerator */
            xmin,               /* minimum pixel x-coord. */
            ymin,               /* minimum pixel y-coord. */
            xmax,               /* maximum pixel x-coord. */
            ymax;               /* maximum pixel y-coord. */

  *min=data[0][0];  *max=data[0][0];
  xmin=0;  ymin=0;
  xmax=0;  ymax=0;
  for (i=0; i<ysize; i++)
  {
    for (j=1; j<xsize; j++)
    {
      if (data[i][j] < *min) { *min=data[i][j]; xmin=j; ymin=i; }
      if (data[i][j] > *max) { *max=data[i][j]; xmax=j; ymax=i; }
    }
  }

  if (vopt)
  {
    printf("min= %lf at (%d, %d)\t\t", *min, xmin, ymin);
    printf("max= %lf at (%d, %d)\n",   *max, xmax, ymax);
  }

  return;
}
/*--------------------------------------------------------*/
int calc_mean(double *dat, int n, double *mean, double *sdev)
{
        int     i;
        double  s,
                ss;

  s=0.0;
  ss=0.0;
  for (i=0; i<n; i++)
  {
    s+=dat[i];
    ss+=dat[i]*dat[i];
  }
  *mean=s/n;
  *sdev=sqrt((ss-s*s/n)/n);

  return(n);
}
/*--------------------------------------------------------*/
int calc_mean2(double minsky, double maxsky, double **dat, int n1, int n2,
               double *mean, double *sdev)
{
        int     i,
                j,
                nn;
        double  s,
                ss;

  s=0.0;
  ss=0.0;
  nn=0;
  for (i=0; i<n2; i++)
  {
    for (j=0; j<n1; j++)
    {
      if ((dat[i][j] > minsky) && (dat[i][j] < maxsky))
      {
        nn++;
        s+=dat[i][j];
        ss+=dat[i][j]*dat[i][j] ;
      }
    }
  }
  *mean=s/nn;
  *sdev=sqrt((ss-s*s/nn)/nn);

  return(nn);
}
/*--------------------------------------------------------*/
void make_hist(double **data, int xsize, int ysize, double min, double max,
      int hs, double bin_width, int *hist_buf)
{
        int i,
            j,
            hi;

  memset((void *)hist_buf, 0, hs*sizeof(int));

  for (i=0; i<ysize; i++)
  {
    for (j=0; j<xsize; j++)
    {
      if ((data[i][j] > min) && (data[i][j] < max))
      {
        hi=(int)((data[i][j]-min)/bin_width);
        if (hi < 0)   p_errmess("make_hist(): hist_buf -> hi < 0");
        if (hi > hs)  p_errmess("make_hist(): hist_buf -> hi > hsize");
        hist_buf[hi]++;
      }
    }
  }

  return;
}
/*--------------------------------------------------------*/
double calc_median(int *hist_buf, int hsize, int limit)
{
        int   i,
              k;
        double median;

  k=0;
  for (i=0; i<hsize; i++)
  {
    k+=hist_buf[i];
    if (k > limit) break;
  }

  median=i-1.0+((double)limit-k+hist_buf[i])/hist_buf[i];

  return(median);
}
/*--------------------------------------------------------*/
double qlsq(size_t n, double *data)
{
        size_t  i;                /* loop numerator         */
        double  s0, s1, s2, s3, s4,
                a, b, c, d, e, f,
                sy, sxy, sxxy,
                xo, yo,
                *nx, *ny,
                det,
                pb,
                pc,
                vertex,
                tmp;

  if (!(nx=(double *)calloc(n, sizeof(double)))) errmess("calloc(nx)");
  if (!(ny=(double *)calloc(n, sizeof(double)))) errmess("calloc(ny)");

  s1 = s2 = 0.0;
  for (i=0; i<n; i++)
  {
    s1+=i;
    s2+=data[i];
  }
  xo=s1/n;
  yo=s2/n;

  for (i=0; i<n; i++)
  {
    nx[i]=(double)(i-xo);
    ny[i]=data[i]-yo;
  }

  s0 = s1 = s2 = s3 = s4 = 0.0;
  pb = pc = 0.0;

  for (i=0; i<n; i++)
  {
    s0+=1.0;

    tmp=nx[i];
    s1+=tmp;

    tmp*=nx[i];
    s2+=tmp;

    tmp*=nx[i];
    s3+=tmp;

    tmp*=nx[i];
    s4+=tmp;
  }

  det=s0*s2*s4+2*s1*s2*s3-s2*s2*s2-s0*s3*s3-s1*s1*s4;

  a=(s2*s4-s3*s3)/det;
  b=(s2*s3-s1*s4)/det;
  c=(s1*s3-s2*s2)/det;
  d=(s0*s4-s2*s2)/det;
  e=(s1*s2-s0*s3)/det;
  f=(s0*s2-s1*s1)/det;

  sy = sxy = sxxy = 0.0;

  for (i=0; i<n; i++)
  {
    sy+=ny[i];
    sxy+=nx[i]*ny[i];
    sxxy+=nx[i]*nx[i]*ny[i];
  }

  pb=b*sy+d*sxy+e*sxxy;
  pc=c*sy+e*sxy+f*sxxy;

  free(nx);
  free(ny);

  vertex=xo-pb/pc/2.0;

  return(vertex);
}
/*--------------------------------------------------------*/
double calc_mode(PARAMS par, size_t hsize, int *hist_buf)
{
        int     hmax;   /*  hist_buf maximum               */
        size_t  i,
                m1,
                m2,
                mbs,
                mi;     /*  hist_buf maximum index         */
        double   mode,
                *mbuf;

  mi=0;
  hmax=hist_buf[0];
  for (i=1; i<hsize; i++)
  {
    if (hist_buf[i] > hmax)
    {
      hmax=hist_buf[i];
      mi=i;
    }
  }
  if (par.vopt > 2)
    printf("calc_mode: histogram maximum at %ld = %d\n", (long)mi, hmax);

  for (i=mi; i>0; i--) if (hist_buf[i] <= hmax/2) break;
  m1=i+1;
  for (i=mi; i<hsize; i++) if (hist_buf[i] <= hmax/2) break;
  m2=i-1;
  mbs=m2-m1+1;
  if (par.vopt > 2)
    printf("calc_mode:  m1= %ld  m2= %ld  mbs= %ld\n",
            (long)m1, (long)m2, (long)mbs);

  if (mbs < 5)
  {
    if (par.vopt) printf("WARNING: Cannot calculate mode.\n");
    return(0.0);
  }
  else
  {
    if (!(mbuf=(double *)calloc(mbs, sizeof(double))))
    {
      fprintf(stderr, "calc_mode: min= %ld  max= %ld => size(mbuf)= %ld\n",
                        (long)m1, (long)m2, (long)mbs);
      errmess("calloc(mbuf)");
    }
    for (i=0; i<mbs; i++) mbuf[i]=(double)hist_buf[i+m1];

    mode=qlsq(mbs, mbuf)+m1;
    if (par.vopt > 2) printf("calc_mode: mode= %lf\n", mode);
    free(mbuf);
  }

  return(mode);
}
/*--------------------------------------------------------*/
double stats(PARAMS par, int xsize, int ysize, double **data)
{
        int     hsize,          /* histogram buffer size  */
                size,           /* number of pixels       */
                nn,             /* num. of pixels for mean */
                *hist_buf;      /* histogram buffer       */
        double   min,            /* minimum pixel value    */
                max,            /* maximum pixel value    */
                bin_width,      /* histogram buffer bin width */
                mode,           /* mode of counts         */
                median,         /* median of counts       */
                skyl,           /* sky level              */
                sigm,           /* std. dev. of counts    */
                mean;           /* mean count rate        */

  size=xsize*ysize;  
  minmax(par.vopt, data, xsize, ysize, &min, &max);

  nn=calc_mean2(par.minsky, par.maxsky, data, xsize, ysize, &mean, &sigm);
  if (par.vopt)  printf("mean=   %G +- %G  (%d pixels)\n", mean, sigm, nn);

  if (par.maxsky > max)
  {
    par.maxsky=max;
    if (par.vopt) printf("maxsky changed to %g\n", max);
  }

  if (par.minsky < min)
  {
    par.minsky=min;
    if (par.vopt) printf("minsky changed to %g\n", min);
  }

  bin_width=0.5;                        /*  bin_width=sigm/BIN_SCALE; */
  hsize=(par.maxsky-par.minsky)/bin_width+1;
  if (!(hist_buf=(int *)calloc(hsize, sizeof(int))))
  {
    printf("sky minimum count= %g   sky maximum count= %g   hist_len= %d\n",
            par.minsky, par.maxsky, hsize);
    errmess("calloc(hist_buf)");
  }
  make_hist(data, xsize, ysize, par.minsky, par.maxsky, hsize, bin_width,
            hist_buf);

  median=calc_median(hist_buf, hsize, size/2)*bin_width+par.minsky;
  if (par.vopt) printf("median= %G\n", median);

  if (fabs(mode=calc_mode(par, hsize, hist_buf)) > 1.0e-10)
    mode=mode*bin_width+par.minsky;
  else
    mode=median;
  if (par.vopt) printf("mode=   %G\n\n", mode);

  free(hist_buf);

  skyl=3.0*median-2.0*mean; /* Stetson */
  if (!strcmp(par.sky_mod, "mode"))
  {
    if (isnan(mode) == 1)
    {
      if (isnan(median) == 1)
      {
        if (isnan(mean) == 1) p_errmess("stats(): cannot find skyl\n");
        else skyl=mean;
      }
      else skyl=median;
    }
    else skyl=mode;
  }
  else if (!strcmp(par.sky_mod, "median"))
  {
    if (isnan(median) == 1)
    {
      if (isnan(mean) == 1)
      {
        if (isnan(mode) == 1) p_errmess("stats(): cannot find skyl\n");
        else skyl=mode;
      }
      else skyl=mean;
    }
    else skyl=median;
  }
  else if (!strcmp(par.sky_mod, "mean"))
  {
    if (isnan(mean) == 1)
    {
      if (isnan(median) == 1)
      {
        if (isnan(mode) == 1) p_errmess("stats(): cannot find skyl\n");
        else skyl=mode;
      }
      else skyl=median;
    }
    else skyl=mean;
  }
  else p_errmess("stats(): unknown algorithm for sky level");

  return(skyl);
}
/*--------------------------------------------------------*/
/** remove multiple points for the same star **/
int reject(short vopt, double minsky, double **data, int n, int *xco, int *yco)
{
        int     *newx,
                *newy,
                *group,
                gn,
                cgn,
                m,
                i,
                j;
        double   *countlim;

  if (!(group=(int *)calloc(n, sizeof(int))))
    errmess("calloc(group)");

  for (i=0; i<n; i++) group[i]=0;

  gn=0;
  for (i=0; i<n; i++)
  {
    cgn=0;
    if (group[i] == 0)
    {
      for (j=0; j<n; j++)
      {
        m=yco[j]-yco[i];
        if (m < -1) continue;
        if (m > 1) break;

        if (j == i) continue;

        if (group[j] != 0)
        {
          if ((abs(xco[i]-xco[j]) < 2) && (abs(yco[i]-yco[j]) < 2))
          {
            cgn=group[j];
            break;
          }
        }
      }

      if (cgn == 0)
      {
        gn++;
        group[i]=gn;
        cgn=gn;
      }
    }
    else cgn=group[i];

    m=0;
    for (j=0; j<n; j++)
    {
      if (j == i) continue;
      if ((abs(xco[i]-xco[j]) < 2) && (abs(yco[i]-yco[j]) < 2))
      {
        m++;
      	if (group[j] == 0) group[j]=cgn;
/*      	else if (group[j]!=cgn)
      	{
      	  printf("i=%d  j=%d  group[i]=%d  group[j]=%d\n",
                  i, j, group[i], group[j]);
          p_errmess("reject(): problem with groups");
        } */
      }
    }
  }

  if (vopt) printf("%d groups of bright pixels\n", gn);

  if (!(newx=(int *)calloc(gn, sizeof(int))))
    errmess("reject: calloc(newx)");
  if (!(newy=(int *)calloc(gn, sizeof(int))))
    errmess("reject: calloc(newy)");
  if (!(countlim=(double *)calloc(gn, sizeof(double))))
    errmess("reject: calloc(countlim)");

  for (i=0; i<gn; i++) countlim[i]=minsky;

  for (i=1; i<=gn; i++)
  {
    for (j=0; j<n; j++)
    {
      if (group[j] == i)
      {
        if (data[yco[j]][xco[j]] > countlim[i-1])
        {
          countlim[i-1]=data[yco[j]][xco[j]];
          newx[i-1]=xco[j];
          newy[i-1]=yco[j];
        }
      }
    }
  }

  for (i=0; i<gn; i++)
  {
    xco[i]=newx[i];
    yco[i]=newy[i];
  }

  free(countlim);
  free(group);
  free(newx);
  free(newy);

  return(gn);
}
/*--------------------------------------------------------*/
int find(PARAMS par, double **data, double clim, int xsize, int ysize,
         int **xco, int **yco)
{
        int i,                  /* loop numerator         */
            j,                  /* loop numerator         */
            *lxco,              /* stars x-coords table   */
            *lyco,              /* stars y-coords table   */
            k,                  /* num. of stars found    */
            n;                  /* num. of bright pixels  */

  if (!(lxco=(int *)calloc(1, sizeof(int)))) errmess("find: calloc(lxco)");
  if (!(lyco=(int *)calloc(1, sizeof(int)))) errmess("find: calloc(lyco)");

  n=0;
  for (i=par.margin; i<ysize-par.margin; i++)
  {
    for (j=par.margin; j<xsize-par.margin; j++)
    {
      if ((data[i][j] > clim) && (data[i][j] < par.maxpeak))
      {
        n++;
        if (!(lxco=(int *)realloc(lxco, n*sizeof(int))))
          errmess("find: realloc(lxco)");
        if (!(lyco=(int *)realloc(lyco, n*sizeof(int))))
          errmess("find: realloc(lyco)");
        lxco[n-1]=j;
        lyco[n-1]=i;
      }
    }
  }

  if (par.vopt)  printf("%d bright pixels\n", n);

  k=reject(par.vopt, par.minsky, data, n, lxco, lyco);

  if (!(lxco=(int *)realloc(lxco, k*sizeof(int))))
    errmess("find: realloc(lxco)");
  if (!(lyco=(int *)realloc(lyco, k*sizeof(int))))
    errmess("find: realloc(lyco)");

  *xco=lxco;
  *yco=lyco;

  if (par.vopt)  printf("%d possible stars\n\n", k);

  return(k);
}
/*--------------------------------------------------------*/
void centroid(PARAMS par, double **data, int x, int y, int xsize, int ysize,
              double skyl, double *cx, double *cy)
{
        int     i,
                j;
        double  tmp,
                sx,
                sy,
                sc;

  if (par.vopt > 2)
    printf("centroid: x= %d/%d  y= %d/%d  aperture= %d\n",
                      x, xsize, y, ysize, par.aperture);

  sx = sy = sc = 0.0;
  for (i=-par.aperture; i<=par.aperture; i++)
  {
    for (j=-par.aperture; j<=par.aperture; j++)
    {
      if (i*i+j*j < par.aperture*par.aperture)
      {
        tmp=data[y+i][x+j]-skyl;
        sx+=j*tmp;
        sy+=i*tmp;
        sc+=tmp;
      }
    }
  }

  *cx=sx/sc+x;
  *cy=sy/sc+y;

  return;
}
/*--------------------------------------------------------*/
char plsq(short vopt, double *x, double *y, double *w, int n, double *a)
{
        int     ai;
        double  den,
                sw,
                sx,
                sy,
                sxx,
                sxy;

  sw = sx = sy = sxx = sxy = 0.0;

  for (ai=0; ai<n; ai++)
  {
    sw+=w[ai];
    sx+=x[ai]*w[ai];
    sy+=y[ai]*w[ai];
    sxx+=x[ai]*x[ai]*w[ai];
    sxy+=x[ai]*y[ai]*w[ai];
  }

  den=sw*sxx-sx*sx;
  if (den == 0.0)
  {
    if (vopt)  printf("Singular matrix.\n");
    return(0);
  }

  a[0]=(sy*sxx-sx*sxy)/den;
  a[1]=(sw*sxy-sy*sx)/den;

  return(1);
}
/*--------------------------------------------------------*/
double find_local_sky(PARAMS par, double **data, int xc, int yc,
                      int xsize, int ysize)
{
        char    mopt;
        int     i, j,
                m,
                hsize,
                *hist_buf;
        double   **lsd,
                bin_width,
                smin, smax,
                r1,
                r2,
                dist,
                sigm,
                sky;

  mopt=par.vopt;
  par.vopt=0;
  if (mopt > 2)
    printf("find_local_sky: x= %d/%d  y= %d/%d\n", xc, xsize, yc, ysize);

  r1=par.iskyr*par.maxfwhm;
  r2=par.oskyr*par.maxfwhm;

/* create local data table */
  if (!(lsd=(double **)calloc(1, sizeof(double *))))
    errmess("calloc(lsd)");

  if (!(lsd[0]=(double *)calloc(1, sizeof(double))))
    errmess("calloc(lsd[0])");

  m=0;
  for (i=-r2; i<=r2; i++)
  {
    for (j=-r2; j<=r2; j++)
    {
      dist=sqrt(i*i+j*j);
      if (dist > r1 && dist < r2)
      {
      	if ((yc+j < 0) || (yc+j > ysize-1)) continue;
      	if ((xc+i < 0) || (xc+i > xsize-1)) continue;
      	
        if (!(lsd[0]=(double *)realloc(lsd[0], (m+1)*sizeof(double))))
          errmess("realloc(lsd[0])");

        lsd[0][m]=data[yc+j][xc+i];
        m++;
      }
    }
  }

  if (mopt > 2) printf("find_local_sky: m= %d\n", m);

  minmax(par.vopt, lsd, m, 1, &smin, &smax);
  if (mopt > 2) printf("find_local_sky: min= %g  max= %g\n", smin, smax);

  bin_width=1.0;                        /*  bin_width=sigm/BIN_SCALE; */
  hsize=(smax-smin)/bin_width+1;
  if (!(hist_buf=(int *)calloc(hsize, sizeof(int))))
  {
    printf("sky minimum count= %g   sky maximum count= %g   hist_len= %d\n",
            smin, smax, hsize);
    errmess("calloc(hist_buf)");
  }
  make_hist(lsd, m, 1, smin, smax, hsize, bin_width, hist_buf);

  if (!strcmp(par.sky_mod, "mode"))
  {
    sky=smin+calc_mode(par, hsize, hist_buf)*bin_width;
  }
  else if (!strcmp(par.sky_mod, "median"))
  {
    sky=calc_median(hist_buf, hsize, m/2)*bin_width+smin;
  }
  else if (!strcmp(par.sky_mod, "mean"))
  {
    calc_mean2(par.minsky, par.maxsky, lsd, m, 1, &sky, &sigm);
  }
  else p_errmess("find_local_sky(): unknown algorithm for sky level");

  par.vopt=mopt;

  free(hist_buf);
  free(lsd[0]);
  free(lsd);

  return(sky);
}
/*--------------------------------------------------------*/
char fitgauss(PARAMS par, double **data, double skyl, int xsize, int ysize, int n,
              double cx, double cy, double *a, int ap)
{
        int     i,
                j,
                k,
                k1;
        double  *w,
                *gx,
                *gy,
                tmp,
                dist;

  if (!(w=(double *)calloc(1, sizeof(double))))
    errmess("calloc(w)");
  if (!(gx=(double *)calloc(1, sizeof(double))))
    errmess("calloc(gx)");
  if (!(gy=(double *)calloc(1, sizeof(double))))
    errmess("calloc(gy)");

  k=0;
  for (i=cy-ap; i<=cy+ap; i++)
  {
    for (j=cx-ap; j<=cx+ap; j++)
    {
      dist=(cx-j)*(cx-j)+(cy-i)*(cy-i);

      if (dist <= 0.0) continue;

      if (sqrt(dist) > ap)  continue;

      tmp=data[i][j]-skyl;
      if (tmp <= 0.0) continue;

      k1=k+1;
      if (!(w=(double *)realloc(w, k1*sizeof(double))))
        errmess("realloc(w)");
      if (!(gx=(double *)realloc(gx, k1*sizeof(double))))
        errmess("realloc(gx)");
      if (!(gy=(double *)realloc(gy, k1*sizeof(double))))
        errmess("realloc(gy)");

      w[k]=sqrt(tmp)/dist;

      gx[k]=dist;
      gy[k]=log(tmp);

      k++;
    }
  }

  if (plsq(par.vopt, gx, gy, w, k, a) == 0)
  {
    free(w);
    free(gx);
    free(gy);

    return(0);
  }

  free(w);
  free(gx);
  free(gy);

  if (a[1] >= 0.0)
  {
    if (par.vopt) printf("Can't fit Gaussian\n");
    return(0);
  }

  a[0]=exp(a[0]);
  a[1]=2.0*sqrt(-log(2.0)/a[1]);

  if (isnan(a[0]) || a[0] < par.minpeak || a[0] > par.maxpeak) return(0);
  if (isnan(a[1]) || a[1] < par.minfwhm || a[1] > par.maxfwhm) return(0);

  return(1);
}
/*--------------------------------------------------------*/
void sort(double *data, int n)
{
        int   i,
              j;
        double tmp;

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
int local_skies(PARAMS par, double **data, int xsize, int ysize, int nstar,
                int **xco, int **yco, double **skyls)
{
        int   i,
              *lxco,
              *lyco,
              n;
        double ts,
              *lskyl;

  if (par.vopt > 1) printf("Determining local sky levels.\n");

  if (!(lskyl=(double *)calloc(nstar, sizeof(double))))
    errmess("calloc(lskyl)");
  lxco=*xco;
  lyco=*yco;

  n=0;
  for (i=0; i<nstar; i++)
  {
    if (par.vopt > 2) printf("local_skies: %d  %d %d\n", i, lxco[i], lyco[i]);
    ts=find_local_sky(par, data, lxco[i], lyco[i], xsize, ysize);
    if (isnan(ts) == 1 || ts < par.minsky || ts > par.maxsky ) continue;

    lskyl[n]=ts;
    lxco[n]=lxco[i];
    lyco[n]=lyco[i];
    n++;
  }

  if (!(lskyl=(double *)realloc(lskyl, n*sizeof(double))))
  {
    printf("Number of sky data points= %d\n", n);
    errmess("realloc(lskyl)");
  }
  if (!(lxco=(int *)realloc(lxco, n*sizeof(int))))
    errmess("realloc(lxco)");
  if (!(lyco=(int *)realloc(lyco, n*sizeof(int))))
    errmess("realloc(lyco)");

  *skyls=lskyl;

  if (par.vopt > 1)
  {
    for (i=0; i<n; i++)
      printf("%5d  (%d,%d) \t %.1f\n", i, lxco[i], lyco[i], lskyl[i]);
    printf("--------------------------------\n");
  }

  return(n);
}
/*--------------------------------------------------------*/
int det_cent(PARAMS par, double **data, int xsize, int ysize, int nstar,
             int *xco, int *yco, double *lskyl, double **cx, double **cy)
{
        int   i,
              j,
              k,
              n;
        double ccx,              /* tmp centroid x-coord.  */
              ccy,              /* tmp centroid y-coord.  */
              *lcx,             /* x-coords of centroids  */
              *lcy,             /* y-coords of centroids  */
              d,
              dx,
              dy;

  if (par.vopt > 1) printf("determining centroids\n");

  if (!(lcx=(double *)calloc(nstar, sizeof(double))))
    errmess("calloc(lcx)");
  if (!(lcy=(double *)calloc(nstar, sizeof(double))))
    errmess("calloc(lcy)");

  if (par.vopt > 2) printf("lcx and lcy allocated, nstar= %d\n", nstar);

  n=0;
  for (i=0; i<nstar; i++)
  {
    if (par.vopt > 3) printf("i= %d/%d\n", i, nstar);

    if ((xco[i] < par.aperture) || (xco[i] > xsize-par.aperture)) continue;
    if ((yco[i] < par.aperture) || (yco[i] > ysize-par.aperture)) continue;

    centroid(par, data, xco[i], yco[i], xsize, ysize, lskyl[i], &ccx, &ccy);
    if (par.vopt > 3) printf("ccx= %lf  ccy= %lf\n", ccx, ccy);
    if ((fabs(ccx-xco[i]) < 1.0) && (fabs(ccy-yco[i]) < 1.0))
    {
      lcx[n]=ccx;
      lcy[n]=ccy;
      n++;
    }
  }

  for (i=0; i<n-1; i++)
  {
    for (j=i+1; j<n; j++)
    {
      dx=lcx[i]-lcx[j];
      dy=lcy[i]-lcy[j];
      d=sqrt(dx*dx+dy*dy);

      if (d < par.minfwhm)
      {
        lcx[i]=(lcx[i]+lcx[j])/2.0;
        lcy[i]=(lcy[i]+lcy[j])/2.0;

        for (k=j; k<n-1; k++)
        {
          lcx[k]=lcx[k+1];
          lcy[k]=lcy[k+1];
        }
        n--;
      }
    }
  }

  if (!(lcx=(double *)realloc(lcx, n*sizeof(double))))
    errmess("realloc(lcx)");
  if (!(lcy=(double *)realloc(lcy, n*sizeof(double))))
    errmess("realloc(lcy)");

  *cx=lcx;
  *cy=lcy;

  if (par.vopt > 1)
  {
    for (i=0; i<n; i++)
      printf("%5d  (%.1f,%.1f)\n", i, lcx[i], lcy[i]);
    printf("--------------------------------\n");
  }
  if (par.vopt) printf("centroids determined for %d stars\n", n);

  return(n);
}
/*--------------------------------------------------------*/
int fwhms(PARAMS par, double **data, int xsize, int ysize, int nstar,
      double *cx, double *cy, double *ffwhm, double *lskyl, double *medf, char *output)
{
        char    *sit;
        int     i,              /* loop numerator         */
                k,
                n,
                nr,
                difap,          /* apertures difference   */
                ap,             /* aperture               */
                oldap;          /* old aperture           */
        double   *fwhm,          /* individual FWHMs table */
                *tsf,           /* FWHM table for stats   */
                *peak,          /* individual PEAKs table */
                a[2],           /* 0 -> PEAK, 1 -> FWHM   */
                mf,             /* mean FWHM              */
                sdev;           /* std.dev. of FWHM       */
	FILE  * outf;


  if (!(sit=(char *)calloc(nstar, sizeof(char))))
    errmess("calloc(sit)");

  if (!(fwhm=(double *)calloc(nstar, sizeof(double))))
    errmess("calloc(fwhm)");
  if (!(tsf=(double *)calloc(nstar, sizeof(double))))
    errmess("calloc(tsf)");
  if (!(peak=(double *)calloc(nstar, sizeof(double))))
    errmess("calloc(peak)");

  *medf=-99.99;
  difap=0;
  ap=par.aperture;
  do
  {
    if (par.vopt) printf("aperture= %d\n", ap);

    oldap=ap;
    for (i=0; i<nstar; i++) sit[i]=1;

    n=0;
    for (i=0; i<nstar; i++)
    {
      if (par.vopt > 2) printf("fitgauss: i= %d  (%lf, %lf)\n", i, cx[i], cy[i]);

      if ((cx[i]< ap) || (xsize-cx[i]< ap) || (cy[i]< ap) || (ysize-cy[i]< ap))
      {
        sit[i]=0;
        fwhm[i]=0.0;
        peak[i]=0.0;
        continue;
      }

      if (!fitgauss(par, data, lskyl[i], xsize, ysize, i, cx[i], cy[i], a, ap))
      {
        sit[i]=0;
        fwhm[i]=0.0;
        peak[i]=0.0;
        continue;
      }

      fwhm[i]=a[1];
      peak[i]=a[0];
      n++;

      if (par.vopt > 1) printf("%5d  (%.1f,%.1f) \t %.1f  %.1f  %.2f\n",
                            i, cx[i], cy[i], lskyl[i], a[0], a[1]);
    }
    if (par.vopt > 1) printf("--------------------------------\n");

    if (n < 2)
    {
      if (par.vopt) printf("Only %d star(s) for statistics\n", n);
      return(-99.99);
    }

    for (k=0; k<9; k++)
    {
      n=0;
      for (i=0; i<nstar; i++)
      {
        if (sit[i])
        {
          tsf[n]=fwhm[i];
          n++;
        }
      }

      calc_mean(tsf, n, &mf, &sdev);

      if (par.vopt > 1) printf("<fwhm>=%lf+-%lf\n", mf, sdev);
      if (sdev < 0.1) break;

      nr=0;
      for (i=0; i<nstar; i++)
      {
        if (!sit[i]) continue;

        if (fabs(fwhm[i]-mf) >= 3.0*sdev)
        {
          sit[i]=0;
          nr++;
          if (par.vopt > 1)
            printf("rejected: %5d  (%.1f,%.1f) \t %.1f  %.2f\n",
                              i, cx[i], cy[i], peak[i], fwhm[i]);
        }
      }
      if (nr == 0) break;
    }

    if (par.vopt > 1) printf("%d average iteration(s)\n", k+1);

    if (n < 2)
    {
      if (par.vopt) printf("Only %d stars in statistics\n", n);
      return(-99.99);
    }

    if (isnan(mf) == 1) return(-99.99);

    ap=(int)mf;
    if (difap*(ap-oldap) < 0) break;
    difap=ap-oldap;
  }
  while(ap!=oldap);

  sort(tsf, n);
  *medf=tsf[n/2];
  

  if (par.vopt)
  {
    for (i=0; i<nstar; i++)
    {
      if (sit[i])
        printf("%3d  ( %6.1f, %6.1f )  sky=%6.1f  peak= %7.1f  FWHM= %5.2f\n",
                  i, cx[i]+1.0, cy[i]+1.0, lskyl[i], peak[i], fwhm[i]);
    }
    printf("--------------------------------\n");
    printf("%d stars  <FWHM>=%lf +- %lf  median=%lf\n\n", n, mf, sdev, *medf);
  }


  if (!(outf = fopen(output, "w"))) errmess("outfname");
  for (i=0;i<nstar;i++)
  {
    if (sit[i])
  	fprintf(outf, "%6.1f %6.1f %6.1f %7.1f %5.2f\n", cx[i]+1.0, cy[i]+1.0, lskyl[i], peak[i], fwhm[i]);
  }
  fclose(outf);

  if (n < 3)  printf("  only %d stars taken ", n);

  *ffwhm=mf;

  free(sit);
  free(tsf);
  free(fwhm);
  free(peak);

  return(n);
}
/*--------------------------------------------------------*/
int main(int argc, char *argv[])
{
        char      *iname,       /* input file name        */
                  *parfname,
                  *ofname,
                  **header;     /* FITS header            */
        int       i,            /* loop numerator         */
                  *xco,         /* stars x-coordinates    */
                  *yco,         /* stars y-coordinates    */
                  hsize,        /* number of header lines */
                  nstar,        /* number of stars        */
                  naxis1,       /* image x-size           */
                  naxis2;       /* image y-size           */
        double     fwhm,         /* mean FWHM              */
                  medf,         /* median FWHM            */
                  skyl,         /* sky counts level       */
                  *lskyl,       /* local sky levels table */
                  *cx,
                  *cy,
                  thresh,       /* finding threshold (above skyl) */
                  tmplim,       /* finding threshild (counts)     */
                  **data;       /* data table             */
        PARAMS    par;

  iname=NULL;
  parfname=NULL;

  if (argc != 4) usage();

  //char * a = argv[1];
  //char * b = argv[2];
  //char * c = argv[3];
  //printf("%s",a);
  //printf("%s",b);
  //printf("%s",c);
  parfname = argv[1];
  iname = argv[2];
  ofname = argv[3];

  printf("%s",parfname);
  printf("%s",iname); 
  printf("%s",ofname); 

  par=readpar(parfname);

  if (!strcmp(par.datatype, "unsigned"))
    data=read_FITS_2Dfile(iname, 'u', &hsize, &header, &naxis1, &naxis2);
  else
    data=read_FITS_2Dfile(iname, 's', &hsize, &header, &naxis1, &naxis2);

  for (i=0; i<hsize; i++) free(header[i]);
  free(header);

  skyl=stats(par, naxis1, naxis2, data);
  if (par.vopt) printf("global sky level= %lf\n\n", skyl);

  nstar=0;
  thresh=par.minpeak;
  for (i=0; i<40; i++)
  {
    tmplim=skyl+thresh;
    if (par.vopt) printf("i=%d  threshold= %lf -> ", i, tmplim);
    nstar=find(par, data, tmplim, naxis1, naxis2, &xco, &yco);

    if (nstar > 100000)  thresh+=par.minpeak;
    else if (nstar < 5)
    {
      if (thresh>par.minpeak) thresh-=par.minpeak/2.0;
      else                    break;
    }
    else break;
  }

  if (par.vopt > 1)
  {
    for (i=0; i<nstar; i++)
      printf("%5d  (%d,%d)\n", i, xco[i], yco[i]);
    printf("--------------------------------\n");
  }

  if (nstar < 2)
  {
    if (par.vopt) printf("Too few stars found\n");
    printf("%s -99.999  -99.999  -99.999  0\n", iname);
    exit(EXIT_FAILURE);
  }

  nstar=local_skies(par, data, naxis1, naxis2, nstar, &xco, &yco, &lskyl);

  nstar=det_cent(par, data, naxis1, naxis2, nstar, xco, yco, lskyl, &cx, &cy);
  free(xco);
  free(yco);

  nstar=fwhms(par, data, naxis1, naxis2, nstar, cx, cy, &fwhm, lskyl, &medf, ofname);

  for (i=0; i<naxis2; i++) free(data[i]);
  free(data);
  free(cx);
  free(cy);
  free(lskyl);

  if (par.vopt) printf("frame       <FWHM>   skylevel  medianFWHM  NSTARS\n");
  printf("%s   %.2f   %.2f   %.2f    %d\n", iname, fwhm, skyl, medf, nstar);

  return(0);
}
/*** END ***/
