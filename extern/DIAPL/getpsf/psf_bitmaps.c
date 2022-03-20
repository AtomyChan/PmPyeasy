/*========================================================*/
/*                                                        */
/*  psf_bitmaps.c       version 1.3.1   2005.01.28        */
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

#include <stdlib.h>

#include "defs.h"
#include "errmess.h"

void psf_bitmaps(IMAGE *im, OBJ_LIST *psf_objects, PSF_STRUCT *psf)
{
        int     npsf, ipsf, cx, cy, buf, size, hw, i, j, nx;
        float   *pix;
        double  *psf_buf;
        OBJECT  *obj;

  nx  = im->nx;
  pix = im->pix;

  obj  = psf_objects->list;
  npsf = psf_objects->nobj;

  hw   = psf->hw;
  size = 2*hw+1;
  buf  = size*size;

  if (!(psf_buf = (double *)malloc(npsf*buf*sizeof(double))))
    errmess("malloc(psf_buf)");
  psf->buf = psf_buf;

  for (ipsf=0; ipsf<npsf; ipsf++)
  {
    cx  = obj->x;
    cy  = obj->y;

    for (i=-hw; i<=hw; i++)
      for (j=-hw; j<=hw; j++)
        psf_buf[hw+i+size*(hw+j)] = pix[cx+i+nx*(cy+j)];

    psf_buf += buf;

    ++obj;
  }

  return;
}
