/*========================================================*/
/*                                                        */
/*  gasdev.c            version 1.2.4   2005.02.16        */
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

#include "ran1.h"

float gasdev(int *idum)
{
        static int    iset=0;
        static float  gset;

        float         fac, r, v1, v2;

  if (iset == 0)
  {
    do
    {
      v1=2.0*ran1(idum)-1.0;
      v2=2.0*ran1(idum)-1.0;
      r=v1*v1+v2*v2;
    }
    while ((r >= 1.0) || (r == 0.0));

    fac=sqrt(-2.0*log(r)/r);
    gset=v1*fac;
    iset=1;
    return(v2*fac);
  }
  else
  {
    iset=0;
    return(gset);
  }
}
