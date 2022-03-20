/*========================================================*/
/*                                                        */
/*  indexx.c            version 1.3.0   2006.04.04        */
/*                                                        */
/*  Modifications:                                        */
/*  Copyright (C) 2006 by Wojtek Pych, CAMK PAN           */
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

void indexx(int n, float *arrin, int *indx)
{
        int   l, j, ir, indxt, i;
        float q;

  for (j=1; j<=n; j++)  indx[j]=j;

  if (n < 2) return;

  l=(n >> 1) + 1;
  ir=n;

  for (;;)
  {
    if (l > 1)  q=arrin[(indxt=indx[--l])];
    else
    {
      q=arrin[(indxt=indx[ir])];
      indx[ir]=indx[1];
      if (--ir == 1)
      {
        indx[1]=indxt;
        return;
      }
    }

    i=l;
    j=l << 1;

    while ((j <= ir) && (j > 0))
    {
      if ((j < ir) && (arrin[indx[j]] < arrin[indx[j+1]])) j++;
      if (q < arrin[indx[j]])
      {
        indx[i]=indx[j];
        j += (i=j);
        if ((i < 1) || (i > (n+1)) || (j < 1) || (j > (n+1))) break;
      }
      else  j=ir+1;
    }
    indx[i]=indxt;
  }

  return;
}
