/*========================================================*/
/*                                                        */
/*  quick_sort.c        version 1.3.1   2005.01.28        */
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

void quick_sort_1(float *, int *, int, int);

/*--------------------------------------------------------*/
void quick_sort(float *list, int *index, int n)
{
        int i;

  for (i=0; i<n; i++)  index[i]=i;

  quick_sort_1(list, index, 0, n-1);

  return;
}
/*--------------------------------------------------------*/
void quick_sort_1(float *list, int *index, int left_end, int right_end)
{
        int   i, j, temp;
        float chosen;

  chosen = list[index[(left_end + right_end)/2]];
  i = left_end-1;
  j = right_end+1;

  for (;;)
  {
    while (list[index[++i]] < chosen);
    {
      while (list[index[--j]] > chosen);

      if (i < j)
      {
        temp=index [j];
        index [j] = index [i];
        index [i] = temp;
      }
      else if (i == j)
      {
        ++i;
        break;
      }
      else break;
    }
  }

  if (left_end < j)  quick_sort_1(list, index, left_end, j);
  if (i < right_end) quick_sort_1(list, index, i, right_end);

  return;
}
