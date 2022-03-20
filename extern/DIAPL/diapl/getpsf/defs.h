/*========================================================*/
/*                                                        */
/*  defs.h              version 1.3.1   2005.01.28        */
/*  (getpsf)                                              */
/*                                                        */
/*  Original source code:                                 */
/*  Copyright (C) 2005 by Przemek Wozniak                 */
/*  wozniak@lanl.gov                                      */
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
  int x;
  int y;
  double xfit;
  double yfit;
  double max;
  double flux;
  double bg;
  double rat;
  double norm;
  double corr;
  int box_index;
} OBJECT;

typedef struct
{
  int nobj;
  int *index;
  OBJECT *list;
} OBJ_LIST;
 
typedef struct
{
  double *pix;
  int nx;
  int ny;
  double sat_level;
  double min_level;
  double gain;
} IMAGE;

typedef struct
{
  int nbox_x;
  int nbox_y;
  int ndeg_spat;
  int ndeg_local;
  int ngauss;
  int npsf_max;
  int nobj_init;
  int min_nbox;
  double min_flux;
  double max_thresh;
  double rat_thresh;
  double nsig_detect;
  double contrast;
  double nsig_rat;  
  int psfhw;
  int maxhw;
  int peakhw;
  int isohw;
  double iso_off;
  double iso_slo;
  double fitrad;
  double aprad;
  double anrad1;
  double anrad2;
  double bkg_frac;
  double nsig_clip;
  int niter_init;
  int niter;
  int verbose;
} FITPAR;

typedef struct
{
  int *list;
  int *index;
  double *buff;
  int irad;
  int nbkg;
  double frac;
} BKG_STRUCT;

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
  double gain;
  double x_orig;
  double y_orig;
} PSF_STRUCT;
