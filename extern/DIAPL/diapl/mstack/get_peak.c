/*========================================================*/
/*                                                        */
/*  get_peak.c          version 1.4.3   2005.02.07        */
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

void get_peak(double *hist, double *bins, int nlo, int nhi,
              double *low, double *high)
{
        int   i, max_i, lo_i, hi_i;
        double max;

  max = -1.0;
  max_i = 0;

  for (i=nlo; i<=nhi; i++)
  {
    if (hist[i] > max)
    {
      max = hist[i];
      max_i = i;
    }
  }

  i = max_i;
  while ((i >= nlo) && (hist[i] > 0.5*max))  i--;
  lo_i = i;

  i = max_i;
  while ((i <= nhi) && (hist[i] > 0.5*max))  i++;
  hi_i = i;

  *low  = bins[lo_i];
  *high = bins[hi_i];

  return;
}
