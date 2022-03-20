/*****************************************************************************/
/* format.h								     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Simple *printf()-like implementation for arbitrary string formatting	     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Copyright (C) 2007; Pal, A. (apal@szofi.elte.hu)			     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*  This library is free software: you can redistribute it and/or modify     */
/*  it under the terms of the GNU General Public License as published by     */
/*  the Free Software Foundation, either version 3 of the License, or	     */
/*  (at your option) any later version.					     */
/*									     */
/*  This program is distributed in the hope that it will be useful,	     */
/*  but WITHOUT ANY WARRANTY; without even the implied warranty of	     */
/*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the	     */
/*  GNU General Public License for more details.			     */
/*									     */
/*  You should have received a copy of the GNU General Public License	     */
/*  along with the program.  If not, see <http://www.gnu.org/licenses/>.     */
/*****************************************************************************/

#ifndef	__FORMAT_H_INCLUDED
#define	__FORMAT_H_INCLUDED	1

/*****************************************************************************/

#define	FORMAT_STRING		0x01
#define	FORMAT_OCT_INTEGER	0x11
#define	FORMAT_DEC_INTEGER	0x12
#define	FORMAT_HEX_INTEGER	0x13
#define	FORMAT_HXC_INTEGER	0x14
#define	FORMAT_INTEGER		FORMAT_DEC_INTEGER

/*****************************************************************************/

/* call: format_check_if_formatted(fmt,"<c1><c2>..."); */
int	format_check_if_formatted(char *format,char *fchars);

/* call: format_replace(fmt,esc,'c1',TYPE1,value1,'c2',TYPE2,value2,...,0); */
char *	format_replace(char *format,int is_escape,...);

/*****************************************************************************/

#endif	/* #ifndef __FORMAT_H_INCLUDED */

/*****************************************************************************/
                                                                        
          
