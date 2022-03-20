/*========================================================*/
/*                                                        */
/*  ran1.c              version 1.2.4   2005.02.16        */
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

#include <stdio.h>
#include <stdlib.h>

#define M1  259200
#define IA1 7141
#define IC1 54773
#define RM1 (1.0/M1)
#define M2  134456
#define IA2 8121
#define IC2 28411
#define RM2 (1.0/M2)
#define M3  243000
#define IA3 4561
#define IC3 51349

double ran1(int *idum)
{
        static int    iff=0;
        static long   ix1, ix2, ix3;
        static double  r[98];

        int   j;
        double temp;

  if ((*idum < 0) || (iff == 0))
  {
    iff=1;
    ix1=(IC1-(*idum)) % M1;
    ix1=(IA1*ix1+IC1) % M1;
    ix2=ix1 % M2;
    ix1=(IA1*ix1+IC1) % M1;
    ix3=ix1 % M3;
    for (j=1; j<=97; j++)
    {
      ix1=(IA1*ix1+IC1) % M1;
      ix2=(IA2*ix2+IC2) % M2;
      r[j]=(ix1+ix2*RM2)*RM1;
    }
    *idum=1;
  }

  ix1=(IA1*ix1+IC1) % M1;
  ix2=(IA2*ix2+IC2) % M2;
  ix3=(IA3*ix3+IC3) % M3;
  j=1 + ((97*ix3)/M3);
  if ((j > 97) || (j < 1))
  {
    fprintf(stderr, "\nRAN1: This cannot happen.\n\n");
    exit(1);
  }

  temp=r[j];
  r[j]=(ix1+ix2*RM2)*RM1;

  return(temp);
}

#undef M1
#undef IA1
#undef IC1
#undef RM1
#undef M2
#undef IA2
#undef IC2
#undef RM2
#undef M3
#undef IA3
#undef IC3
