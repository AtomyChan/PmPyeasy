/*========================================================*/
/*                                                        */
/*  hedit.c             version 1.3.0   2005.03.03        */
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
/*                                                        */
/* Edits keywords in FITS header header[hsize][80].       */
/* name is the name of keyword (full 8 characters),       */
/* type=c, s, i, d, f for character, string, int, double  */
/*      and double values,                                 */
/* descr is comment to appear after keyword value         */
/*       and / sign,                                      */
/* fmt is format to be used (just like in printf),        */
/* and the last argument (...) is the value of keyword.   */
/* In case name is not found in the header,               */
/* or type of keyword is not recognized, the header       */
/* remains unchanged.                                     */
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
#include <stdarg.h>
#include <string.h>

void hedit(int hsize, char **header, char *name, char type, char *descr,
           char *fmt, ...)
{
        char    cval, *sval, entry[256];
        int     len, i, ival, found;
        double  dval;
        va_list ap;

  va_start(ap, fmt);

  switch(type)
  {
    case 'c': cval = (char)va_arg(ap, int);
              len  = sprintf(entry, fmt, cval);
              break;
    case 's': sval = va_arg(ap, char *);
              len  = sprintf(entry, fmt, sval);
              break;
    case 'i': ival = va_arg(ap, int);
              len  = sprintf(entry, fmt, ival);
              break;
    case 'd':
    case 'f': dval = va_arg(ap, double);
              len  = sprintf(entry, fmt, dval);
              break;
    default:  fprintf(stderr, "hedit: warning: unknown type of keyword\n");
              fprintf(stderr, "header unchanged\n");
              len = 0;
              break;
  }

  if (len != 0)
  {
    len += sprintf(entry+len, " / %s", descr);

    for (i=len; i<70; i++) entry[i] = ' ';

    found = 0;
    for (i=0; i<hsize; i++)
    {
      if (strncmp(header[i], name, 8) == 0)
      {
        memcpy(header[i]+10, entry, 70);
        found = 1;
      }
    }

    if (!found)
    {
      fprintf(stderr, "hedit: warning: %s not found\n", name);
      fprintf(stderr, "header unchanged\n");
    }
  }

  return;
}
