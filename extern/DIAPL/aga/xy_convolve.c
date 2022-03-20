/*========================================================*/
/*                                                        */
/*  xy_convolve.c       version 1.3.0   2005.03.03        */
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

#include <malloc.h>

void xy_convolve(float *im, double *cnvim, double *tmp, int nx, int ny,
                 double *ker_x, double *ker_y, int kerhw)
{
        int i, j, k, idx;

  for (j=0; j<ny; j++)
  {
    for (i=kerhw; i<nx-kerhw; i++)
    {
      idx = i + nx*j;
      tmp[idx] = 0.0;

      for (k=-kerhw; k<=kerhw; k++)
        tmp[idx] += im[k+idx]*ker_x[kerhw-k];
    }
  }

  for (i=kerhw; i<nx-kerhw; i++)
  {
    for (j=kerhw; j<ny-kerhw; j++)
    {
      idx = i + nx*j;
      cnvim[idx] = 0.0;

      for (k=-kerhw; k<=kerhw; k++)
        cnvim[idx] += tmp[idx+nx*k]*ker_y[kerhw-k];
    }
  }

  return;
}
