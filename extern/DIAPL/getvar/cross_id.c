/*========================================================*/
/*                                                        */
/*  cross_id.c          version 1.2.3   2005.01.24        */
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

#include "defs.h"

void cross_id(STAR *obj1, int nobj1, STAR *obj2, int nobj2, float id_rad)
{
        int   i, j;
        float rr, rad2, dx, dy;

  rad2 = id_rad*id_rad;

  for (j=0; j<nobj2; j++) obj2[j].vtype = 10;

  for (i=0; i<nobj1; i++)
  {
    obj1[i].vtype = 1;

    for (j=0; j<nobj2; j++)
    {
      dx = obj2[j].x - obj1[i].x;
      dy = obj2[j].y - obj1[i].y;

      rr = dx*dx + dy*dy;

      if (rr < rad2)
      {
        obj2[j].x = 0.0;
        obj2[j].y = 0.0;
        obj2[j].nframes = 0;
        obj1[i].vtype += 10;
        break;
      }
    }
  }

  return;
}
