/*========================================================*/
/*                                                        */
/*  tri_gen.c           version 1.3.4   2006.04.03        */
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
/* Triangle generation. From the list of nobj points(x,y) */
/* calculate perimeters, longest/shortest side ratios,    */
/* cosines of the angle between longest/shortest sides    */
/* and orientations for all possible triangles.           */
/* Returns also indexes of vertices ver1-2-3 opposite     */
/* to the shortest, mid and longest side and number       */
/* of triangles.                                          */
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

#include "indexx.h"

void tri_gen(float *x, float *y, int nobj, float *perim, float *rat, float *cos,
             int *ori, int *ver1, int *ver2, int *ver3, long *ntri)
{
        int   i, j, k, idx[3], vi[3];
        long  n;
        float xc, yc, dx, dy, r[3], s[3];

  n=0;

  for (i=0; i<nobj; i++)
  {
    for (j=i+1; j<nobj; j++)
    {
      for (k=j+1; k<nobj; k++)
      {
	vi[0]=k;
	vi[1]=j;
	vi[2]=i;

	dx=x[j]-x[i];
	dy=y[j]-y[i];
	r[0]=sqrt(dx*dx+dy*dy);

	dx=x[k]-x[i];
	dy=y[k]-y[i];
	r[1]=sqrt(dx*dx+dy*dy);

	dx=x[k]-x[j];
	dy=y[k]-y[j];
	r[2]=sqrt(dx*dx+dy*dy);

	xc=(x[i]+x[j]+x[k])/3;
	yc=(y[i]+y[j]+y[k])/3;

/* use little trick to fool stupid 1,n indexing of numerical recipes */
	indexx(3, r-1, idx-1);

	ver1[n]=vi[idx[0]-1];
	ver2[n]=vi[idx[1]-1];
	ver3[n]=vi[idx[2]-1];

	s[0]=r[idx[0]-1];
	s[1]=r[idx[1]-1];
	s[2]=r[idx[2]-1];

	ori[n]=(x[ver2[n]]-xc)*(y[ver1[n]]-yc);
	ori[n]-=(x[ver1[n]]-xc)*(y[ver2[n]]-yc);
	ori[n]/=fabs(ori[n]);

	perim[n]=s[2];
	rat[n]=s[2]/s[0];
	cos[n]=(s[0]*s[0]+s[2]*s[2]-s[1]*s[1])/(2.0*s[0]*s[2]);

	n++;
      }
    }
  }

  *ntri=n;

  return;
}
