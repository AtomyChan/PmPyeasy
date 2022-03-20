/*========================================================*/
/*                                                        */
/*  read_sector.c       version 1.4.0   2005.04.19        */
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

#include "errmess.h"
#include "pfitsio1.h"

void read_sector(char *infn, long headlen, int nx0, int ny0,
                int x0_off, int y0_off, char *im, int bpx, int nx, int ny)
{
        int   j, bpr;
        long  ofs;
        FILE  *inf;

  if (!(inf = fopen(infn, "rb"))) errmess(infn);

  ofs = headlen + (x0_off + nx0*y0_off)*bpx;
  fseek(inf, ofs, SEEK_SET);

  bpr = nx*bpx;
  ofs = nx0*bpx - bpr;

  for (j=0; j<ny; j++)
  {
    fread(im+bpr*j, 1, bpr, inf);
    fseek(inf, ofs, SEEK_CUR);
  }

  if (needswap()) swap_bytes(im, nx*ny, bpx);

  fclose(inf);

  return;
}
