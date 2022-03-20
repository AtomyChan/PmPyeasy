/*========================================================*/
/*                                                        */
/*  median.c            version 1.4.3   2005.02.07        */
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

float median(float *hist, float *bins, int nbin,
             float low, float high, float frac)
{
        int     i, lo_i;
        float   value;
        double  sum1, sum2;

  sum1 = sum2 = 0.0;

  i = 0;
  while ((bins[i] < low) && (i < nbin)) i++;

  lo_i = i;
  while ((bins[i] <= high) && (i < nbin)) sum1 += hist[i++];

  i = lo_i;
  while ((bins[i] <= high) && (sum2 < frac*sum1) && (i < nbin))
    sum2 += hist[i++];

  value = bins[i];

  return(value);
}
