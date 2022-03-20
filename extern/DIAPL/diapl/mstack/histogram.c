/*========================================================*/
/*                                                        */
/*  histogram.c         version 1.4.3   2005.02.07        */
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

void histogram(float *im, int nx, int ny, float *hist, int nbin,
               float low, float high)
{
        int   i, idx;
        float bin_size;

  bin_size = (high - low)/(float)nbin;

  for (i=0; i<nbin; i++)  hist[i] = 0.0;

  for (i=0; i<nx*ny; i++)
  {
    idx = (int)((im[i] - low)/bin_size);

    if (idx >= nbin)  idx = nbin-1;
    if (idx < 0)      idx = 0;

    hist[idx] += 1.0;
  }

  return;
}
