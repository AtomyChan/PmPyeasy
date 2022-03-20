/*========================================================*/
/*                                                        */
/*  .c               version 1.3.1   2005.01.28        */
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

void copy_obj(OBJECT *obj1, OBJECT *obj2)
{
  obj1->x         = obj2->x;
  obj1->y         = obj2->y;
  obj1->xfit      = obj2->xfit;
  obj1->yfit      = obj2->yfit;
  obj1->max       = obj2->max;
  obj1->flux      = obj2->flux;
  obj1->bg        = obj2->bg;
  obj1->rat       = obj2->rat;
  obj1->norm      = obj2->norm;
  obj1->box_index = obj2->box_index;

  return;
}
