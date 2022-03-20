/*========================================================*/
/*                                                        */
/*  bkg.c               version 1.3.1   2005.01.28        */
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

void init_bkg(IMAGE *im, BKG_STRUCT *bkg, float rad1, float rad2)
{
        int   i, j, nbkg, irad;
        float rad, r1, r2;

  r1 = rad1*rad1;
  r2 = rad2*rad2;

  nbkg = 0;
  irad = (int)ceil(rad2);

  if (!( bkg->list = (int *)malloc((2*irad+1)*(2*irad+1)*sizeof(int))))
    errmess("malloc(bkg->list)");

  for (i=-irad; i<=irad; i++)
  {
    for (j=-irad; j<=irad; j++)
    {
      rad = i*i + j*j;

      if ((rad > r1) && (rad <= r2))  bkg->list[nbkg++] = i + im->nx*j;
    }
  }

  if (!(bkg->buff = (float *)malloc(nbkg*sizeof(float))))
    errmess("malloc(bkg->buff)");
  if (!(bkg->index = (int *)malloc(nbkg*sizeof(int))))
    errmess("malloc(bkg->index)");

  bkg->nbkg = nbkg;
  bkg->irad = irad;

  return;
}
/*--------------------------------------------------------*/
double getbkg(IMAGE *im, BKG_STRUCT *bkg, int cx, int cy)
{
        int     k, ofs, irad, nbkg;
        double  bg;

  bg   = 0.0;
  ofs  = cx + im->nx*cy;
  irad = bkg->irad;
  nbkg = bkg->nbkg;

  if ((cx >= irad) && (cx < im->nx-irad) && (cy >= irad) && (cy < im->ny-irad))
  {
    for (k=0; k<nbkg; k++)  bkg->buff[k] = im->pix[bkg->list[k]+ofs];

    quick_sort(bkg->buff, bkg->index, nbkg);

    bg = bkg->buff[bkg->index[(int)(nbkg*bkg->frac)]];
  }

  return(bg);
}
