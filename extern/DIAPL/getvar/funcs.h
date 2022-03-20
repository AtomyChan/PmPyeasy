/*========================================================*/
/*                                                        */
/*  funcs.h             version 1.3.1   2005.04.20        */
/*  (getvar)                                              */
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

extern void   get_phot(float *, float *, float *, double, double, int, int,
                        double *, STAR *, PAR_STRUCT *);
extern void   get_params(char *, char *, PAR_STRUCT *);
extern void   read_psf(char *, PSF_STRUCT *, int);
extern void   init_psf(PSF_STRUCT *, double, double);
extern double psf_core(PSF_STRUCT *, double, double);
extern void   read_kernel(char *, KER_STRUCT *, int, int);
extern void   spatial_coeffs(KER_STRUCT *, double, double, double *);
extern void   make_vectors(KER_STRUCT *);
extern void   make_kernel(KER_STRUCT *, double *, double *, int);
extern void   im_convolve(double *, double *, int, int, double *, int);
extern float  *rfits_f(char *, char *, int *, int *,
                        float *, float *, long int *);
extern void   make_psf(PSF_STRUCT *, double, double, int, int, double *);
extern void   indexx(int, float *, int *);
extern void   base_func(float, int, int, double *, double *, int, int);
extern double get_fwhm(double *, double, double, int, int, PAR_STRUCT *,
                        double *);
extern void   covar_sig(float *, float *, int, int, float, int, double *, int);
extern float  bkg(float *, int, int, PAR_STRUCT *);
extern float  aperphot(float *, int, int, int, int, float, float, int *, float);
extern void   centroid(float *, int, int, int, int, float *, float *);
extern STAR   *find_stars(float *, int, int, int *, int *, PAR_STRUCT *);
extern void   cross_id(STAR *, int, STAR *, int, float);
extern void   quick_sort(float *, int *, int);
extern void   get_repeaters(float **, float **, float *, float *,
                            PAR_STRUCT *, int);
extern void   center_stars(STAR *, int, float **, float **, PAR_STRUCT *, int);
extern int    neighbor(float *, int, int, PAR_STRUCT *, float);


