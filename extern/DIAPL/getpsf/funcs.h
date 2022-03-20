/*========================================================*/
/*                                                        */
/*  funcs.h             version 1.3.2   2005.04.07        */
/*  (getpsf)                                              */
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

extern void get_params(char *, char *, IMAGE *, FITPAR *, PSF_STRUCT *);
extern void init_bkg(IMAGE *, BKG_STRUCT *, float, float);
extern double getbkg(IMAGE *, BKG_STRUCT *, int, int);
extern int is_peak(IMAGE *, int, int, int, float);
extern void local_max(IMAGE *, FITPAR, OBJ_LIST *);
extern void psf_stars(IMAGE *, OBJ_LIST *, OBJ_LIST *, FITPAR);
extern float *rfits_f(char *, char *, int *, int *,
		      float *, float *, long int *);
extern double silly_phot(IMAGE *, OBJECT *, int);
extern void swap(char *, long int, int);
extern int needswap(void);
extern void nrerror(char *);
extern float *vector(int, int);
extern void free_vector(float *, int, int);
extern void refine_max(IMAGE *, OBJECT *);
extern void center(IMAGE *, OBJECT *);
extern int is_isolated(IMAGE *, OBJECT *, FITPAR);
extern void quick_sort(float *, int *, int);
extern void copy_obj(OBJECT *, OBJECT *);
extern void psf_bitmaps(IMAGE *, OBJ_LIST *, PSF_STRUCT *);
extern void fitpsf(OBJ_LIST *, PSF_STRUCT *);
extern void ludcmp(double **, int, int *, double *);
extern void lubksb(double **, int, int *, double *);
extern void get_gaussian(PSF_STRUCT *);
extern void init_psf(PSF_STRUCT *, double, double);
extern double psf_core(PSF_STRUCT *, double, double);
extern void local_fit(OBJ_LIST *, PSF_STRUCT *);
extern void clip_stars(OBJ_LIST *, float, int *);
extern void write_psf(int, char *, PSF_STRUCT *);
/*** END ***/
