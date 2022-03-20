/*========================================================*/
/*                                                        */
/*  funcs.h             version 1.4.2  2005.07.08         */
/*  (aga)                                                 */
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

#include <stdio.h>

/* preparation of arrays for image fitting and clipping */

extern void get_params(char *, char *, PARAMS *par);
extern void expand_matrix(double **, double *, double **, DOM_TABS *, PARAMS);
extern void make_vectors(float *, double **, double **, double **, PARAMS);
extern void make_domains(float *, unsigned short *, double **, DOM_TABS *, PARAMS);
extern void clone_domains(DOM_TABS *, int, int);
extern void mask_badpix(float *, unsigned short *, unsigned short *, int, PARAMS);
extern int  clean_domains(float *, float *, unsigned short *, unsigned short *,
			  double **, DOM_TABS *, PARAMS);
extern void get_domains(float *, unsigned short *, DOM_TABS *, PARAMS, int *);
extern int  max(int, int, int, float *, int);

/* sigma clipping and local value of fit */

extern void quick_sort(float *, int *, int);
extern  int local_clip(float *, float *, unsigned short *, double **,
		       double **, double *, DOM_TABS *, PARAMS, double *, int *);
extern  int clip_domains(DOM_TABS *, float *, int *, int, double*, PARAMS);
extern void local_solution(int, int, double **, double *, double *, PARAMS);

/* basis function for kernel decomposition */

extern double base_func(float, int, int, double *, double *, int, int);

/* FITS manipulation */

extern void init_difimages(char **, char **, long int *, int, PARAMS par);

extern void hedit(char *, int, char *, char, char *, char *, ...);
extern long get_headlen(char *);
extern void read_sector(char *infn, long, int, int, int, int, char *,
                         int, int, int);
extern int  write_sector(char *infn, long int, int, int, int, int, char *,
                         int, int, int, int);

/* convolutions */

extern void xy_convolve(float *, double *, double *, int, int,
			double *, double *, int);
extern void spatial_convolve(float *, float *, double **, double **,
			     double *, PARAMS);

/* solving of linear equation (numerical recipes) */

extern void ludcmp(double **, int, int *, double *);
extern void lubksb(double **, int, int *, double *);

/* output */

extern void write_kernel(char *, double *, int, int, int, PARAMS, double *);
extern void apply_norm(float *, unsigned short *, unsigned short *,
		       double, int, int, float);

