/*========================================================*/
/*                                                        */
/*  defs.h              version 1.2.3   2005.01.24        */
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

typedef struct
{
  int nx0;
  int ny0;
  int psfhw;
  int psfn;
  int bkg_mode;
  int verbose;
  double aprad;
  double fitrad;
  double normrad;
  double anrad1;
  double anrad2;
  double nsig_bkg;
  double sat_level;
  double min_level;
  double bad_value;
  double gain;
  double err;
  double c_min;
  int smhw;
  int mohw;
  int nobj_init;
  int ncons_var1;
  int npts_var2;
  double nsig_var1;
  double nsig_var2;
  double lim_ratio;
  double id_rad;
  int isohw;
  double iso_off;
  double iso_slo;
  double bad_margin;
  int center_hw;
  int filter_hw;
  double fwhm_frac;
} PAR_STRUCT;

typedef struct
{
  double x;
  double y;
  int cx;
  int cy;
  double flux;
  int nframes;
  double *bmp;
  double p_flux;
  double a_flux;
  double p_err;
  double a_err;
  double corr;
  double chi2_n;
  double bg;
  int nbad;
  int vtype;
} STAR;


typedef struct
{
  int hw;
  int ndeg_spat;
  int ndeg_local;
  int ngauss;
  int recenter;
  double *buf;
  double cos;
  double sin;
  double ax;
  double ay;
  double sigma_inc;
  double sigma_mscale;
  double *vec;
  double *local_vec;
  double fitrad;
  double normrad;
  double gain;
  double x_orig;
  double y_orig;
} PSF_STRUCT;

typedef struct
{
  int nsec_x;
  int nsec_y;
  int nx;
  int ny;
  int hw;
  int wdeg;
  int ncomp;
  int nvecs;
  int ntot;
  int nwxy;
  double *sig;
  int *deg;
  double *chi2_n;
  double **vec;
  double **vecs;
} KER_STRUCT;

/*** END ***/
