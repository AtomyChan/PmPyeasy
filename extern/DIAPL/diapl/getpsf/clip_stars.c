/*========================================================*/
/*                                                        */
/*  clip_stars.c        version 1.3.1   2005.01.28        */
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

#define BIG_SIGMA 100.0

void clip_stars(OBJ_LIST *psf_objects, double nsigma, int *nstars)
{
        int     npsf, ipsf, ngood, nclip;
        double  mean, sigma, corr;
        OBJECT  *obj;

  npsf = psf_objects->nobj;
  obj  = psf_objects->list;

  sigma = mean = 0.0;
  ngood = 0;

  for (ipsf=0; ipsf<npsf; ipsf++)
  {
    corr = obj->corr;

    if (corr > 0.0)
    {
      mean  += corr;
      sigma += corr*corr;
      ++ngood;
    }
    ++obj;
  }

  if (ngood > 0)
  {
    mean /= ngood;
    if (ngood > 1)  sigma = sqrt(sigma/ngood - mean*mean);
    else            sigma = BIG_SIGMA;
  }
  else
  {
    sigma = BIG_SIGMA;
    mean  = BIG_SIGMA;
  }

  obj  = psf_objects->list;

  ngood = nclip = 0;

  for (ipsf=0; ipsf<npsf; ipsf++)
  {
    corr = obj->corr;

    if (mean - corr > nsigma*sigma)
    {
      obj->corr = -1.0;
      ++nclip;
    }
    else  ++ngood;

    ++obj;
  }

  *nstars = ngood;

  return;
}
