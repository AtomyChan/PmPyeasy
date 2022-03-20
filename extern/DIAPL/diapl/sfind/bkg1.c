/*========================================================*/
/*                                                        */
/*  bkg.c               version 1.3.3   2006.10.23        */
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
/* Given image im(nx,ny) and centroid (cx,cy) estimates   */
/* background taking data from annulus (r1,r2).           */
/* Sigma clipping with n_sig until nothing changes.       */
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

#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <math.h>

#include "defs.h"
#include "errmess.h"
#include "indexx.h"

/*--------------------------------------------------------*/
void minmax(short verbose, double **data, int xsize, int ysize,
            double *min, double *max)
{
        int i,                  /* rows loop numerator    */
            j,                  /* columns loop numerator */
            xmin,               /* minimum pixel x-coord. */
            ymin,               /* minimum pixel y-coord. */
            xmax,               /* maximum pixel x-coord. */
            ymax;               /* maximum pixel y-coord. */

  *min = *max = data[0][0];
  xmin = ymin = xmax = ymax = 0;

  for (i=0; i<ysize; i++)
  {
    for (j=1; j<xsize; j++)
    {
      if (data[i][j] < *min) { *min=data[i][j]; xmin=j; ymin=i; }
      if (data[i][j] > *max) { *max=data[i][j]; xmax=j; ymax=i; }
    }
  }

  if (verbose > 1)
  {
    printf("\t min= %g at (%d, %d)\t", *min, xmin, ymin);
    printf("max= %g at (%d, %d)\n",   *max, xmax, ymax);
  }

  return;
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

  s = ss = 0.0;
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

  bzero((void *)hist_buf, hs*sizeof(int));

  for (i=0; i<ysize; i++)
  {
    for (j=0; j<xsize; j++)
    {
      if ((data[i][j] > min) && (data[i][j] < max))
      {
        hi=(int)((data[i][j]-min)/bin_width);
        if ((hi < 0) || (hi > hs))
        {
          printf("ERROR! make_hist() hi= %d - outside hist_buf[]\n", hi);
          exit(EXIT_FAILURE);
        }
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
double qlsq(double *data, int n)
{
        int   i;                /* loop numerator         */
        double s0, s1, s2, s3, s4,
              a, b, c, d, e, f,
              sy, sxy, sxxy,
              xo, yo,
              *nx, *ny,
              det,
              pb,
              pc,
              vertex,
              tmp;

  if (!(nx=(double *)calloc((size_t)n, sizeof(double)))) errmess("calloc(nx)");
  if (!(ny=(double *)calloc((size_t)n, sizeof(double)))) errmess("calloc(ny)");

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
    nx[i]=i-xo;
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
double calc_mode(PARAMS par, int *hist_buf, int hsize)
{
        int     i,
                m1,
                m2,
                mbs,
                mi,     /*  hist_buf maximum index         */
                hmax;   /*  hist_buf maximum               */
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
  if (par.verbose > 2)
    printf("calc_mode()\n\t histogram maximum: hmax= %d at mi= %d \n",
      hmax, mi);

  for (i=mi; i>0; i--) if (hist_buf[i] <= hmax/2) break;
  m1=i+1;
  for (i=mi; i<hsize; i++) if (hist_buf[i] <= hmax/2) break;
  m2=i-1;
  if (par.verbose > 2) printf("\t m1= %d  m2= %d\n", m1, m2);

  mbs=m2-m1+1;

  if (mbs < 7)
  {
    if (par.verbose)
      printf("WARNING! calc_mode() mbs= %d cannot calculate mode\n", mbs);
      
    return(sqrt(-1.0)); // intentional error
  }
  else
  {
    if (!(mbuf=(double *)calloc((size_t)mbs, sizeof(double))))
    {
      printf("calc_mode() min=%d  max=%d => size(mbuf)=%d\n", m1, m2, mbs);
      errmess("calloc(mbuf)");
    }
    for (i=0; i<mbs; i++) mbuf[i]=(double)hist_buf[i+m1];

    mode=qlsq(mbuf, mbs)+m1;
    if (par.verbose > 2) printf("\t mode= %lf\n", mode);
    free(mbuf);
  }

  return(mode);
}
/*--------------------------------------------------------*/
double bkg_level(int xsize, int ysize, double **data, PARAMS par, double *sigma)
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
                mean,           /* mean count rate        */
                lminlev, lsatlev;

  size=xsize*ysize;
  lminlev=par.min_level;
  lsatlev=par.sat_level;
  minmax(par.verbose, data, xsize, ysize, &min, &max);

  nn=calc_mean2(lminlev, lsatlev, data, xsize, ysize, &mean, sigma);
  if (par.verbose > 1)
    printf("\t mean=   %G +- %G  (%d pixels)\n", mean, *sigma, nn);

  if (lsatlev > max) lsatlev=max;
  if (lminlev < min) lminlev=min;

  bin_width=1.0;                        /*  bin_width=sigma/BIN_SCALE; */
  hsize=(int)((lsatlev-lminlev)/bin_width)+1;
  if (par.verbose > 2) printf("\t hsize= %d\n", hsize);
  if (!(hist_buf=(int *)calloc((size_t)hsize, sizeof(int))))
    errmess("calloc(hist_buf)");
  make_hist(data, xsize, ysize, lminlev, lsatlev, hsize, bin_width,
            hist_buf);

  median=calc_median(hist_buf, hsize, size/2)*bin_width+lminlev;
  if (par.verbose > 1) printf("\t median= %G\n", median);

  mode=lminlev+calc_mode(par, hist_buf, hsize)*bin_width;
  if (par.verbose > 1) printf("\t mode=   %G\n\n", mode);

  free(hist_buf);

  skyl=3.0*median-2.0*mean; /* Stetson */
  if (par.bkg_algorithm == 1)       // mode
  {
    if (isnan(mode))
    {
      if (isnan(median))
      {
        if (isnan(mean))
        {
          printf("ERROR! stats() cannot find background level\n");
          exit(EXIT_FAILURE);
        }
        else skyl=mean;
      }
      else skyl=median;
    }
    else skyl=mode;
  }
  else if (par.bkg_algorithm == 2)  // median
  {
    if (isnan(median))
    {
      if (isnan(mean))
      {
        if (isnan(mode))
        {
          printf("ERROR! stats() cannot find background level\n");
          exit(EXIT_FAILURE);
        }
        else skyl=mode;
      }
      else skyl=mean;
    }
    else skyl=median;
  }
  else if (par.bkg_algorithm == 3)  // mean
  {
    if (isnan(mean))
    {
      if (isnan(median))
      {
        if (isnan(mode))
        {
          printf("ERROR! stats() cannot find background level\n");
          exit(EXIT_FAILURE);
        }
        else skyl=mode;
      }
      else skyl=median;
    }
    else skyl=mean;
  }
  else
  {
    printf("ERROR! stats() unknown algorithm for background\n");
    exit(EXIT_FAILURE);
  }

  return(skyl);
}
/*--------------------------------------------------------*/
double bkg(int nx, int ny, double **im, int cx, int cy, double r1, double r2,
          PARAMS par)
{
        int     i, j,
                ir,
                *idx,
                delta, lx, ly;
        size_t  npix;
        double   rr, rr1, rr2,
                bglev,
                absdev,
                mean,
                *tab,
                pix;

  npix = (size_t)(M_PI*((r2+1.0)*(r2+1.0) - r1*r1));
  if (par.verbose > 2)
    printf("bkg()\n\t x= %d  y= %d  r1= %g  r2= %g  npix= %ld\n",
            cx, cy, r1, r2, (long)npix);

  if (!(tab=(double *)calloc(npix, sizeof(double))))  errmess("calloc(tab)");
  if (!(idx=(int *)calloc(npix, sizeof(int))))       errmess("calloc(idx)");

  ir=(int)r2+1;
  rr1=r1*r1;
  rr2=r2*r2;

/* load pixels within annulus (r1, r2) into tab */

  npix=0;

  for (j=-ir; j<=ir; j++)
  {
    for (i=-ir; i<=ir; i++)
    {
      rr=(double)i*i+(double)j*j;

      if ((rr1 < rr) && (rr < rr2))
      {
      	lx=cx+i;
      	ly=cy+j;
      	if ((lx >= 0) && (lx < nx) && (ly >= 0) && (ly < ny))
      	{
      	  pix=im[ly][lx];

          if ((pix >= par.min_level) && (pix < par.sat_level))
            tab[npix++]=pix;
        }
      }
    }
  }

/* take median as first background level                             */
/* use little trick to fool stupid 1,n indexing of numerical recipes */

  if (npix < 3) return(0.0);

  indexx((int)npix, tab-1, idx-1);

  bglev=tab[idx[npix/2]-1];

/* reject pixels deviating by more than n_sig*sigma       */
/* from current bglev                                     */
/* subsequent bglev will be means of the remaining pixels */
/* quit when no more rejections                           */

  delta=1;

  while (delta > 0)
  {
    absdev=0.0;

    for (i=0; i<(int)npix; i++) absdev+=(bglev - tab[i])*(bglev - tab[i]);

    absdev*=par.n_sig*par.n_sig/(double)npix;

    j=0;
    mean=0.0;

    for (i=0; i<(int)npix; i++)
    {
      if ((bglev - tab[i])*(bglev - tab[i]) < absdev)
      {
        tab[j++]=tab[i];
        mean+=tab[i];
      }
    }

    if (j > 1)
    {
      delta=(int)npix-j;
      npix=(size_t)j-1;
      mean/=npix;

      bglev=mean;
    }
    else
      return(0.0);
  }

  free(tab);
  free(idx);

  return(bglev);
}
