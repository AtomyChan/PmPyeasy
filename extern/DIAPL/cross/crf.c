/*========================================================*/
/*                                                        */
/*  crf.c               version 1.3.3   2005.01.28        */
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
#include <stdio.h>

#include "defs.h"
#include "funcs.h"

float crf(float *data1, float *data2, float *imout, int *nn, int flag)
{
        int   i, j;
        float Re, Im, sum_Re, sum_Im, rat;
       
  fourn(data2, nn, 2, 1);

  for (i=0; i<2*nn[0]*nn[1]; i+=2)
  {
    j = i + 1;
    imout[i] = data1[i]*data2[i] + data1[j]*data2[j];
    imout[j] = data1[j]*data2[i] - data1[i]*data2[j];
  }

  fourn(imout, nn, 2, -1);

  rat = 0.0;

  if (flag)
  {
    sum_Re = sum_Im = 0.0;

    for (i=0; i<2*nn[0]*nn[1]; i+=2)
    {
      j = i + 1;
      Re = imout[i];
      Im = imout[j];
      sum_Re += Re*Re;
      sum_Im += Im*Im;
    }

    rat = sqrt(sum_Re/sum_Im);
  }

  return(rat);
}

