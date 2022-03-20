/*========================================================*/
/*                                                        */
/*  mask_badpix.c       version 1.3.0   2005.03.03        */
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

#include "defs.h"

void mask_badpix(float *im, unsigned short *mask,
		 unsigned short *flag, int mode, PARAMS par)
{
        int i, j, k, l, idx1, idx2, xlo, xhi, ylo, yhi, bad_growrad;

  if (mode)
  {
    for (idx1=0; idx1<par.nx*par.ny; idx1++)
    {
      mask[idx1] = 1;
      flag[idx1] = 1;
    }
    bad_growrad = par.bad_growrad2;
  }
  else
  {
    for (idx1=0; idx1<par.nx*par.ny; idx1++)
      flag[idx1] = 1;
    bad_growrad = par.bad_growrad1;
  }

  for (i=0; i<par.nx; i++)
  {
    for (j=0; j<par.ny; j++)
    {
      idx1 = i + par.nx*j;

      if ((im[idx1] < par.min_level) || (im[idx1] >= par.sat_level))
      {
	xlo = (i > bad_growrad) ? i-bad_growrad : 0;
	ylo = (j > bad_growrad) ? j-bad_growrad : 0;

	xhi = (i < par.nx-bad_growrad) ? i+bad_growrad : par.nx-1;
	yhi = (j < par.ny-bad_growrad) ? j+bad_growrad : par.ny-1;
	
       	for (k=xlo; k<=xhi; k++)
        {
	  for (l=ylo; l<=yhi; l++)
          {
	    idx2 = k + par.nx*l;
	    mask[idx2] = 0;
	    flag[idx2] = 0;
	  }
	}
      }
    }
  }

  return;
}
/*** END ***/
