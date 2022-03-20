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
  float aprad;
  float fitrad;
  float normrad;
  float anrad1;
  float anrad2;
  float nsig_bkg;
  float sat_level;
  float min_level;
  float bad_value;
  float gain;
  float err;
  float c_min;
  int smhw;
  int mohw;
  int nobj_init;
  int ncons_var1;
  int npts_var2;
  float nsig_var1;
  float nsig_var2;
  float lim_ratio;
  float id_rad;
  int isohw;
  float iso_off;
  float iso_slo;
  float bad_margin;
  int center_hw;
  int filter_hw;
  float fwhm_frac;
} PAR_STRUCT;

typedef struct
{
  float x;
  float y;
  int cx;
  int cy;
  float flux;
  int nframes;
  float *bmp;
  float p_flux;
  float a_flux;
  float p_err;
  float a_err;
  float corr;
  float chi2_n;
  float bg;
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
  float *sig;
  int *deg;
  double *chi2_n;
  double **vec;
  double **vecs;
} KER_STRUCT;

/*** END ***/
