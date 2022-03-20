/*========================================================*/
/*                                                        */
/*  local_fit.c         version 1.3.1   2005.01.28        */
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

#include <math.h>

#include "defs.h"
#include "funcs.h"

void local_fit(OBJ_LIST *psf_objects, PSF_STRUCT *psf)
{
        int     i, j, psfhw, cx, cy, npsf, ipsf,
                size, cxpsf, cypsf, indx[4];
        double  dx, dy, bg, x_off, y_off,
                *psf_buf, rr, parity, rad2, xpsf, ypsf, ampl, x, y,
                sum, sumpsf, x2_psf, y2_psf, xy_psf,
                pix, pix_psf, v[3], vec[4], *mat[4], m1[4], m2[4], m3[4],
                a, b, c, corr;
        OBJECT  *obj;

  mat[1] = m1;
  mat[2] = m2;
  mat[3] = m3;

  psfhw   = psf->hw;
  psf_buf = psf->buf;
  rad2    = psf->fitrad*psf->fitrad;
  size    = 2*psfhw+1;

  npsf = psf_objects->nobj;
  obj  = psf_objects->list;

/* main loop over psf stars */

  for (ipsf=0; ipsf<npsf; ipsf++)
  {
    for (i=0; i<3; i++)
    {
      vec[i+1] = 0.0;
      for (j=0; j<=i; j++)  mat[i+1][j+1] = 0.0;
    }

    xpsf  = obj->xfit;
    ypsf  = obj->yfit;
    cxpsf = obj->x;
    cypsf = obj->y;

    dx = xpsf - (double)cxpsf;
    dy = ypsf - (double)cypsf;

    bg = obj->bg;

/* get local coefficients of the psf */
    init_psf(psf, xpsf, ypsf);

/* --------------------------------- */

    ampl = sum = sumpsf = x2_psf = y2_psf = xy_psf = 0.0;

    a = b = c = 0.0;

    for (cx=-psfhw; cx<=psfhw; cx++)
    {
      for (cy=-psfhw; cy<=psfhw; cy++)
      {
        x = (double)cx - dx;
        y = (double)cy - dy;

        rr = x*x + y*y;

        if (rr <= rad2)
        {
          pix     = psf_buf[cx+psfhw+size*(cy+psfhw)];
          pix_psf = psf_core(psf,x,y);

          ampl   += pix_psf*(pix-bg)/pix;
          sum    += pix_psf*pix_psf/pix;
          sumpsf += pix_psf;

          a += pix_psf*(pix-bg);
          b += pix_psf*pix_psf;
          c += (pix-bg)*(pix-bg);

          x2_psf += pix_psf*x*x;
          y2_psf += pix_psf*y*y;
          xy_psf += pix_psf*x*y;

          v[0] = pix_psf;
          v[1] = pix_psf*x;
          v[2] = pix_psf*y;

          for (i=0; i<3; i++)
          {
            vec[i+1] += v[i]*(pix-bg)/pix;
            for (j=0; j<=i; j++)  mat[i+1][j+1] += v[i]*v[j]/pix;
          }
        }
      }
    }

    for (i=0; i<3; i++)
      for (j=0; j<=i; j++)
        mat[j+1][i+1] = mat[i+1][j+1];

    ludcmp(mat, 3, indx, &parity);
    lubksb(mat, 3, indx, vec);

/* smart linear recentering */

    if (sum != 0.0)
    {
      ampl = ampl/sum;

      if (psf->recenter && (sumpsf != 0.0))
      {
        x_off = (x2_psf*vec[2] + xy_psf*vec[3])/(fabs(vec[1])*sumpsf);
        y_off = (xy_psf*vec[2] + y2_psf*vec[3])/(fabs(vec[1])*sumpsf);

        if (fabs(x_off) > 1.0)  x_off = 0.0;
        if (fabs(y_off) > 1.0)  y_off = 0.0;

        xpsf += x_off;
        ypsf += y_off;
      }
    }

/* correlation coefficient for goodness of fit */

    if ((b != 0.0) && (c != 0.0)) corr = a/sqrt(b*c);
    else                          corr = -1.0;

    obj->xfit = xpsf;
    obj->yfit = ypsf;

    obj->norm = ampl;
    obj->corr = corr;

    ++obj;
    psf_buf += size*size;
  }

  return;
}
