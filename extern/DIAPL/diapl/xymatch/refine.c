/*========================================================*/
/*                                                        */
/*  refine.c            version 1.3.3   2005.02.10        */
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
/* Given first order transformation in coeffx and coeffy  */
/* transform list2 of nobj2 objects (x2,y2)               */
/* into approximate frame of reference of list1           */
/* of nobj1 objects (x1,y1) and look for close neighbors  */
/* to make a refined matching of lists.                   */
/* TOL is size of matching box and upon return index      */
/* contains subscripts of objects from list2              */
/* indentified in list1.                                  */
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

int refine(double *coeffx, double *coeffy, double *x1, double *y1,
            double *x2, double *y2, int nobj1, int nobj2, int *index, double TOL)
{
        int   i, j, n;
        double xt, yt;

  for (i=0; i<nobj1; i++) index[i]=-1;

  n=0;
  for (j=0; j<nobj2; j++)
  {
    xt = coeffx[0]*x2[j] + coeffx[1]*y2[j] + coeffx[2];
    yt = coeffy[0]*x2[j] + coeffy[1]*y2[j] + coeffy[2];

    for (i=0; i<nobj1; i++)
    {
      if (fabs(xt - x1[i]) < TOL  &&  fabs(yt - y1[i]) < TOL)
      {
        index[i]=j;
        n++;
      }
    }
  }

  return(n);
}
/*** END ***/
