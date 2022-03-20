/*========================================================*/
/*                                                        */
/*  get_headlen.c       version 1.4.0   2005.04.19        */
/*                                                        */
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
#include <string.h>

#include "errmess.h"
#include "pfitsio1.h"

long get_headlen(char *inpfname)
{
        char  head[RECORD_SIZE],
              eoh;
        int   i;
        long  ofs;
        FILE  *inpf;

  if (!(inpf=fopen(inpfname, "r"))) errmess(inpfname);

  eoh=0;
  ofs=0L;
  while (!eoh)
  {
    if (fread(head, RECORD_SIZE, 1, inpf) != 1 )
      errmess("get_headlen: reading header");

    for (i=0; i<RECORD_SIZE; i+=CARD_SIZE)
      if (strncmp(head+i, "END     ", 8) == 0)  eoh=1;

    ofs+=RECORD_SIZE;
  }

  fclose(inpf);

  return(ofs);
}
