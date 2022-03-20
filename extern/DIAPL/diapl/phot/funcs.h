/*========================================================*/
/*                                                        */
/*  funcs.h             version 1.8.0   2006.01.02        */
/*  (phot)                                                */
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

extern void get_phot(double *, double *, double *, double, double, int, int,
		     double *, STAR *, PAR_STRUCT *);
extern void get_params(char *, char *, PAR_STRUCT *);
extern void read_psf(char *, PSF_STRUCT *);
extern void init_psf(PSF_STRUCT *, double, double);
extern double psf_core(PSF_STRUCT *, double, double);
extern void read_kernel(char *, KER_STRUCT *, int, int);
extern void spatial_coeffs(KER_STRUCT *, double, double, double *);
extern void make_vectors(KER_STRUCT *);
extern void make_kernel(KER_STRUCT *, double *, double *, int);
extern void im_convolve(double *, double *, int, int, double *, int);
extern double *rfits_f(char *, char *, int *, int *,
		      double *, double *, long int *);
extern void make_psf(PSF_STRUCT *, double, double, int, int, double *);
extern double bkg(double *, int, int, PAR_STRUCT *);
extern void indexx(int, double *, int *);
extern int get_sector(double *, double *, double *, double *, KER_STRUCT *);
extern void base_func(double, int, int, double *, double *, int, int);
extern double get_fwhm(double *, double, double, int, int,
		       PAR_STRUCT *, double *);
