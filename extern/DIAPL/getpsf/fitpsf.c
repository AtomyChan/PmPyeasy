/*========================================================*/
/*                                                        */
/*  fitpsf.c            version 1.3.1   2005.01.28        */
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
#include <math.h>
#include <stdlib.h>

#include "defs.h"
#include "funcs.h"
#include "errmess.h"

void fitpsf(OBJ_LIST *psf_objects, PSF_STRUCT *psf)
{
        int     i, j, k, l, m, n, psfhw, ipsf, cx, cy, igauss, icomp, ispat,
                npsf, ldeg, sdeg, ncomp, nspat, ngauss, size, cxpsf, cypsf,
                *indx, ntot;
        double  xpsf, ypsf, a1, a2, dx, dy, *v, *wxy, bg, norm, weight, rad2,
                **mat, *vec, **mat0, *vec0, *psf_buf, f, rr, x1, y1, x2, y2,
                parity, sigma2_inc, sin_psf, cos_psf, ax_psf, ay_psf, pix,
                x_orig, y_orig;
        OBJECT  *obj;

  sigma2_inc = psf->sigma_inc*psf->sigma_inc;

  obj        = psf_objects->list;
  npsf       = psf_objects->nobj;
  psfhw      = psf->hw;
  psf_buf    = psf->buf;
  sdeg       = psf->ndeg_spat;
  ldeg       = psf->ndeg_local;
  cos_psf    = psf->cos;
  sin_psf    = psf->sin;
   ax_psf    = psf->ax;
   ay_psf    = psf->ay;
  x_orig     = psf->x_orig;
  y_orig     = psf->y_orig;

  ngauss     = psf->ngauss;
  ncomp      = ngauss*(ldeg+1)*(ldeg+2)/2;
  nspat      = (sdeg+1)*(sdeg+2)/2;
  ntot       = ncomp*nspat;
  size       = 2*psfhw+1;
  rad2       = psf->fitrad*psf->fitrad;

/* get memory */

  if (!(wxy = (double *)malloc(nspat*sizeof(double)))) errmess("malloc(wxy)");

  if (!(v    = (double *)malloc(ncomp*sizeof(double)))) errmess("malloc(v)");
  if (!(vec0 = (double *)malloc(ncomp*sizeof(double)))) errmess("malloc(vec0)");
  if (!(mat0 = (double **)malloc(ncomp*sizeof(double *))))
    errmess("malloc(mat0)");

  for (i=0; i<ncomp; i++)
    if (!(mat0[i] = (double *)malloc(ncomp*sizeof(double))))
      errmess("malloc(mat0[i])");

  if (!(vec  = (double *)malloc(ntot*sizeof(double)))) errmess("malloc(vec)");
  if (!(indx = (int *)malloc(ntot*sizeof(int)))) errmess("malloc(indx)");
  if (!(mat = (double **)malloc(ntot*sizeof(double *)))) errmess("malloc(mat)");

  for (i=0; i<ntot; i++)
  {
    if (!(mat[i] = (double *)malloc(ntot*sizeof(double))))
      errmess("malloc(mat[i])");
    mat[i]--;
  }
  mat--;

  for (i=0; i<ntot; i++)
  {
    vec[i] = 0.0;
    for(j=0; j<ntot; j++) mat[i+1][j+1] = 0.0;
  }

/* main loop over all psf stars */

  for (ipsf=0; ipsf<npsf; ipsf++)
  {
    if (obj->corr < 0.0)
    {
      ++obj;
      psf_buf += size*size;
      continue;
    }

    xpsf  = obj->xfit;
    ypsf  = obj->yfit;
    cxpsf = obj->x;
    cypsf = obj->y;

    dx = xpsf - (double)cxpsf;
    dy = ypsf - (double)cypsf;

    bg   = obj->bg;
    norm = obj->norm;

/********************************/
/***        Local Fit         ***/
/********************************/

    for (i=0; i<ncomp; i++)
    {
      vec0[i] = 0.0;
      for (j=0; j<ncomp; j++) mat0[i][j] = 0.0;
    }

    for (cx=-psfhw; cx<=psfhw; cx++)
    {
      for (cy=-psfhw; cy<=psfhw; cy++)
      {
        pix = psf_buf[cx+psfhw+size*(cy+psfhw)];
        weight = pix/(norm*norm);
        pix -= bg;

        x1 = (double)cx - dx;
        y1 = (double)cy - dy;

        x2 = cos_psf*x1 - sin_psf*y1;
        y2 = sin_psf*x1 + cos_psf*y1;

        rr = ax_psf*x2*x2 + ay_psf*y2*y2;

        if (rr <= rad2)
        {
          icomp = 0;
          for (igauss=0; igauss<ngauss; igauss++)
          {
            f = exp(rr);
            a1 = 1.0;
            for (m=0; m<=ldeg; m++)
            {
              a2 = 1.0;
              for (n=0; n<=ldeg-m; n++)
              {
                v[icomp++] = f*a1*a2;
                a2 *= y1;
              }
              a1 *= x1;
            }
            rr *= sigma2_inc;
          }

          if (icomp != ncomp)
            fprintf(stderr, "fitpsf: warning: wrong number of components!\n");

/***********************/

          for (i=0; i<ncomp; i++)
          {
            vec0[i] += v[i]*pix/norm/weight;
            for (j=0; j<=i; j++)  mat0[i][j] += v[i]*v[j]/weight;
          }
        }
      }
    }

    for (i=0; i<ncomp; i++)
      for (j=0; j<i; j++)
        mat0[j][i] = mat0[i][j];

/********************************/
/***       Global Fit         ***/
/********************************/

    ispat = 0;
    a1 = 1.0;
    for (m=0; m<=sdeg; m++)
    {
      a2 = 1.0;
      for (n=0; n<=sdeg-m; n++)
      {
        wxy[ispat++] = a1*a2;
        a2 *= ypsf - y_orig;
      }
      a1 *= xpsf - x_orig;
    }

    if (ispat != nspat)
      fprintf(stderr, "fitpsf: warning: wrong number of spatial components!\n");

    m = 0;
    for (i=0; i<nspat; i++)
    {
      for (j=0; j<ncomp; j++)
      {
        vec[m] += wxy[i]*vec0[j];

        n = 0;
        for (k=0; k<nspat && n<=m; k++)
        {
          for (l=0; l<ncomp && n<=m; l++)
          {
            mat[m+1][n+1] += wxy[i]*wxy[k]*mat0[j][l];
            n++;
          }
        }
        m++;
      }
    }

    psf_buf += size*size;
    ++obj;
  }

  for (i=1; i<=ntot; i++)
    for (j=1; j<i; j++)
      mat[j][i] = mat[i][j];

  ludcmp(mat, ntot, indx-1, &parity);
  lubksb(mat, ntot, indx-1, vec-1);

  for(i=0; i<ntot; i++) psf->vec[i] = vec[i];

/* done with memory ! */

  free(v);
  free(vec);
  free(vec0);
  free(wxy);
  free(indx);

  mat++;
  for (i=0; i<ntot; i++)  free(++mat[i]);
  free(mat);

  for(i=0; i<ncomp; i++)  free(mat0[i]);
  free(mat0);

  return;
}
