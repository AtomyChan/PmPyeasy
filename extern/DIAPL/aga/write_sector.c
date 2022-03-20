/*========================================================*/
/*                                                        */
/*  write_sector.c      version 1.3.0   2005.03.03        */
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

int write_sector(char *outfn, long int headlen, int nx0, int ny0, int x0_off,
                 int y0_off, char *im, int bpx, int nx, int ny, int marg)
{
        char  *pim;
        int   j, bpr, btw;
        long  ofs;
        FILE  *outf;

  if (!(outf = fopen(outfn, "rb+"))) errmess(outfn);

  ofs = headlen + (x0_off + marg + nx0*(y0_off + marg))*bpx;
  fseek(outf, ofs, SEEK_SET);

  bpr = nx*bpx;
  btw = (nx - 2*marg)*bpx;
  ofs = nx0*bpx - btw;
  pim = im + marg*bpx;

  if (needswap()) swap_bytes(im, nx*ny, bpx);

  for (j=marg; j<ny-marg; j++)
  {
    fwrite(pim+bpr*j, 1, btw, outf);
    fseek(outf, ofs, SEEK_CUR);
  }

  if (needswap()) swap_bytes(im, nx*ny, bpx);

  fclose(outf);

  return(1);
}
