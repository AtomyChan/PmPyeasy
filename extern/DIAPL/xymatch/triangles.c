/*========================================================*/
/*                                                        */
/*  triangles.c         version 1.3.4   2006.04.03        */
/*                                                        */
/*  Original source code:                                 */
/*  Copyright (C) 2006 by Przemek Wozniak                 */
/*  wozniak@lanl.gov                                      */
/*                                                        */
/*  Modifications:                                        */
/*  Copyright (C) 2006 by Wojtek Pych, CAMK PAN           */
/*  pych@camk.edu.pl                                      */
/*                                                        */
/*  Written for GNU project C and C++ Compiler            */
/*                                                        */
/* Given two coordinate lists (x1, y1) and (x2, y2)       */
/* (each nobj long) return index xx of objects            */
/* from list2 identified in list1.                        */
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
#include <stdlib.h>
#include <math.h>

#include "tri_gen.h"
#include "errmess.h"
#include "defs.h"

void triangles(float *x1, float *y1, float *x2, float *y2,
               int nobj1, int nobj2, int *xx, PARAMS par)
{
        int   *ori1, *ori2, *ver11, *ver12, *ver13, *ver21, *ver22, *ver23,
              *vote, *vidx, tver=0, count, found, best,
              *index,
              idx;
        long  ntri1, ntri2, i, j, k;
        float *perim1, *perim2, *cos1, *cos2, *rat1, *rat2,
              dperim, drat, dcos;

  ntri1 = nobj1*(nobj1-1)*(nobj1-2)/6;
  ntri2 = nobj2*(nobj2-1)*(nobj2-2)/6;
  if (par.verbose > 2)
    printf("triangles()\n\t ntri1= %ld  ntri2= %ld\n", ntri1, ntri2);

  if (!(perim1=(float *)malloc(ntri1*sizeof(float)))) errmess("malloc(perim1)");
  if (!(perim2=(float *)malloc(ntri2*sizeof(float)))) errmess("malloc(perim2)");
  if (!(rat1  =(float *)malloc(ntri1*sizeof(float)))) errmess("malloc(rat1)");
  if (!(rat2  =(float *)malloc(ntri2*sizeof(float)))) errmess("malloc(rat2)");
  if (!(cos1  =(float *)malloc(ntri1*sizeof(float)))) errmess("malloc(cos1)");
  if (!(cos2  =(float *)malloc(ntri2*sizeof(float)))) errmess("malloc(cos2)");

  if (!(ori1 =(int *)malloc(ntri1*sizeof(int)))) errmess("malloc(ori1)");
  if (!(ori2 =(int *)malloc(ntri2*sizeof(int)))) errmess("malloc(ori2)");
  if (!(ver11=(int *)malloc(ntri1*sizeof(int)))) errmess("malloc(ver11)");
  if (!(ver12=(int *)malloc(ntri1*sizeof(int)))) errmess("malloc(ver12)");
  if (!(ver13=(int *)malloc(ntri1*sizeof(int)))) errmess("malloc(ver13)");
  if (!(ver21=(int *)malloc(ntri2*sizeof(int)))) errmess("malloc(ver21)");
  if (!(ver22=(int *)malloc(ntri2*sizeof(int)))) errmess("malloc(ver22)");
  if (!(ver23=(int *)malloc(ntri2*sizeof(int)))) errmess("malloc(ver23)");

  if (!(index=(int *)malloc(ntri1*sizeof(int)))) errmess("malloc(index)");

/* generate triangles from both lists of points */
  tri_gen(x1, y1, nobj1, perim1, rat1, cos1, ori1, ver11, ver12, ver13, &ntri1);
  tri_gen(x2, y2, nobj2, perim2, rat2, cos2, ori2, ver21, ver22, ver23, &ntri2);
  if (par.verbose > 2)
    printf("\t ntri1= %ld  ntri2= %ld\n", ntri1, ntri2);

/* for triangles big and "round" enough look for close neighbors */
/* in (perimeter, longest/shortest side ratio, cosine angle      */
/* between longest/shortest side) space                          */
  for (i=0L; i<ntri1; i++)
  {
    index[i]=-1;

    for (j=0L; j<ntri2; j++)
    {
      if ((rat1[i] < par.rlim) && (rat2[j] < par.rlim) &&
          (perim1[i] > par.llim) && (perim2[j] > par.llim))
      {
        dperim=fabs(perim1[i] - perim2[j]);
        drat  =fabs(1.0 - rat1[i]/rat2[j]);
        dcos  =fabs(1.0 - cos1[i]/cos2[j]);

        if ((dperim < par.ltol) && (drat < par.rtol) &&
            (dcos < par.ctol) && (ori1[i] == ori2[j]))
          index[i]=j;
      }
    }
  }

  free(perim1); free(rat1); free(cos1); free(ori1);
  free(perim2); free(rat2); free(cos2); free(ori2);

/* voting procedure: for each point check how many triangles */
/* were matched and had this point as one of its vertices    */

  if (!(vote=(int *)calloc(ntri1, sizeof(int)))) errmess("calloc(vote)");
  if (!(vidx=(int *)calloc(ntri1, sizeof(int)))) errmess("calloc(vidx)");

  for (i=0L; i<nobj1; i++)
  {
    count=0;

    for (j=0L; j<ntri1; j++)
    {
      if ((idx=index[j]) != -1)
      {
        if (i == ver11[j])      tver = ver21[idx];
        else if (i == ver12[j]) tver = ver22[idx];
        else if (i == ver13[j]) tver = ver23[idx];
        else                    idx = -1;
        if (idx != -1)
        {
          found=0;
          for (k=0L; k<count; k++)
          {
            if (vidx[k] == tver)
            {
              vote[k]++;
              found=1;
            }
          }

          if (!found)
          {
            vidx[count]=tver;
            vote[count]=1;
            count++;
          }
        }
      }
    }

    best=0;
    for (k=0L; k<count; k++) if (vote[k] > vote[best]) best=k;
    if (par.verbose > 2)
      printf("%ld count= %d  vote= %d  best= %d  vote(best)= %d (%g)\n",
        i, count, vote[count], best, vote[best], par.fvno*nobj1*(nobj1-1)/2);

    if ((count) && (vote[best] > par.fvno*nobj1*(nobj1-1)/2))
      xx[i]=vidx[best];
    else
      xx[i]=-1;
  }

  free(ver11); free(ver12); free(ver13);
  free(ver21); free(ver22); free(ver23);

  free(vote);
  free(vidx);

  return;
}
/*** END ***/
