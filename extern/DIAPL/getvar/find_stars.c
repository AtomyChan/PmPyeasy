/*========================================================*/
/*                                                        */
/*  find_stars.c        version 1.2.3   2005.01.24        */
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

STAR *find_stars(float *corrim, int nx, int ny, int *nobj, int *nobj_max,
                 PAR_STRUCT *par)
{
        int   i, j, k, l, idx[4], max_test, SMHW, MOHW, LIM, nsig_bkg;
        float lims[4], cpix, sat_level, anrad1, anrad2, aprad,
              gain, c_min;
        STAR  *obj;

  if (!(obj=(STAR *)malloc(par->nobj_init*sizeof(STAR))))
    errmess("malloc(obj)");

  *nobj_max = par->nobj_init;
  *nobj = 0;

  sat_level = par->sat_level;
  nsig_bkg  = (int)(par->nsig_bkg);

  anrad1 = par->anrad1;
  anrad2 = par->anrad2;
  aprad  = par->aprad;
  gain   = par->gain;

  c_min  = par->c_min;
  SMHW   = par->smhw;
  MOHW   = par->mohw;

  lims[0] = par->aprad;
  lims[1] = par->anrad2;
  lims[2] = SMHW;
  lims[3] = 2*MOHW+1;

  indexx(4, lims-1, idx-1);

  LIM = lims[idx[3]-1];

  for (i=LIM; i<nx-LIM; i++)
  {
    for (j=LIM; j<ny-LIM; j++)
    {
      cpix = corrim[i+nx*j];
      max_test = 1;

/* look for local extrema in covariance matrix */

      if (fabs(cpix) >= c_min)
      {
        for (l=-SMHW; l<=SMHW; l++)
          for (k=-SMHW; k<=SMHW; k++)
            if (cpix < corrim[i+k+nx*(j+l)])
              max_test = 0;

/* for local extrema find cetroid, background and photometry */

        if (max_test)
        {
          if ((*nobj) >= (*nobj_max))
          {
            (*nobj_max) += (par->nobj_init);
            if (!(obj = (STAR *)realloc(obj,(*nobj_max)*sizeof(STAR))))
              errmess("find_stars: realloc(obj)");
          }

          obj[*nobj].cx = i;
          obj[*nobj].cy = j;

          (*nobj)++;
        }
      }
    }
  }

  return (obj);
}
