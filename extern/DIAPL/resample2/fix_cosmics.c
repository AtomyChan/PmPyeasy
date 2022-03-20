/*========================================================*/
/*                                                        */
/*  fix_cosmics.c       version 1.2.4   2005.02.16        */
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
#include <stdlib.h>

#include "indexx.h"
#include "errmess.h"

void fix_cosmics(float *im, int nx, int ny,
                 float MAX_GRAD, float N_SIG, float GAIN)
{
           int            i, j, k, l, nmed, *idx, indx;
           float          mean, rat,
                          med[9], med0, medim;

  if (!(idx=(int *)malloc(nx*ny*sizeof(int)))) errmess("malloc(idx)");

  indexx(nx*ny, im-1, idx-1);
  med0 = im[idx[nx*ny/2]-1];

  for (j=1; j<ny-1; j++)
  {
    for(i=1; i<nx-1; i++)
    {
      indx = j*nx + i;
      nmed = 0;
      for (l=-1; l<=1; l++)
        for (k=-1; k<=1; k++)
          med[nmed++] = im[indx+l*nx+k];

      indexx(nmed, med-1, idx-1);
      medim = med[idx[nmed/2]-1];

      mean = medim - med0;

      if (mean < 1.0) mean = 1.0;

      rat = (im[indx] - med0)/mean;

      if ((rat > MAX_GRAD) && (im[indx] > N_SIG*sqrt(im[indx]/GAIN)))
        im[indx] = medim;
    }
  }

  free(idx);

  return;
}
