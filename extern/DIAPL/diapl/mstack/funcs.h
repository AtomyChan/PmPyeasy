/*========================================================*/
/*                                                        */
/*  funcs.h             version 1.6.0   2005.05.06        */
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

extern int read_sector(char *infn, size_t, int, int, int, int,
		       char *, int, int, int);
extern int write_sector(char *infn, long int, int, int, int, int,
			char *, int, int, int, int);
extern void get_params(char *, char *, PAR_STRUCT *);
extern void quick_sort(double *, int *, int);
extern double median(double *, double *, int, double, double, double);
extern void histogram(double *, int, int, double *, int, double, double);
extern void bin(double *, int, double *, int);
extern void get_peak(double *, double *, int, int, double *, double *);
