/*========================================================*/
/*                                                        */
/*  im_convolve.c       version 1.6.2   2005.01.24        */
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

void im_convolve(double *im, double *cnvim, int nx, int ny,
                 double *kerim, int kerhw)
{
        int i, j, m, n, kern, idx;

  kern = 2*kerhw + 1;

  for (j=kerhw; j<ny-kerhw; j++)
  {
    for (i=kerhw; i<nx-kerhw; i++)
    {
      idx = i + nx*j;
      cnvim[idx] = 0.0;

      for (n=-kerhw; n<=kerhw; n++)
        for(m=-kerhw; m<=kerhw; m++)
          cnvim[idx] += im[m+idx+nx*n]*kerim[kerhw-m+kern*(kerhw-n)];
    }
  }

  return;
}
