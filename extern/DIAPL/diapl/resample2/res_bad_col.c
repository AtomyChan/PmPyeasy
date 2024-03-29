/*========================================================*/
/*                                                        */
/*  res_bad_col.c       version 1.2.5   2005.04.07        */
/*                                                        */
/*  Modifications:                                        */
/*  Copyright (C) 2005 by Wojtek Pych, CAMK PAN           */
/*  pych@camk.edu.pl                                      */
/*                                                        */
/*  Written for GNU project C and C++ Compiler            */
/*                                                        */
/* Restore bad column. Function interpolates every        */
/* bad columns onto the grid of a reference image         */
/* and change values of bad pixels to 0.                  */
/* (C) Igor Soszynski & Karol Zebrun                      */
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

#include "poly.h"
#include "defs.h"

void res_bad_col(RINGPAR par, double *im, int nx, int ny, int *col,
                 double *coeffx, int ndeg)
{
        int   i, j,
              npix,
              k;
        double x;

  npix=nx*ny;

  for (i=0; i<nx; i++)
  {
    if (col[i] > par.sat)
    {
      if (par.verbose > 3) printf("res_bad_col: restoring column: %d\n", i);
      for (j=0; j<ny; j++)
      {
        x = poly((double)i, (double)j, ndeg, coeffx);
        if (par.verbose > 3) printf("res_bad_col: x= %g\n", x);

        k=j*nx+(int)(2*i-x); 
        if (par.verbose > 3) printf("res_bad_col: k1= %d\n", k);
        if ((k >= 0) && (k < npix)) im[k] = 0.0;

        k=j*nx+(int)(2*i-x+1); 
        if (par.verbose > 3) printf("res_bad_col: k2= %d\n", k);
        if ((k >= 0) && (k < npix)) im[k] = 0.0;
      }
    }
  }

  return;
}
