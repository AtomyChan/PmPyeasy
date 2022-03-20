/*========================================================*/
/*                                                        */
/*  defs.h              version 1.7.0   2005.05.05        */
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

typedef struct
{
  char  dbf;
  int   nx0, ny0,
        psfhw,
        psfn,
        bkg_mode,
        verbose;
  float aprad,
        fitrad,
        normrad,
        anrad1, anrad2,
        nsig_bkg,
        sat_level,
        min_level,
        bad_value,
        gain,
        err;
} PAR_STRUCT;

typedef struct
{
  int   nbad;
  float p_flux,
        a_flux,
        p_err,
        a_err,
        corr,
        chi2_n,
        ker_chi2,
        bg,
        fwhm;
} STAR;

typedef struct
{
  int     hw,
          ndeg_spat,
          ndeg_local,
          ngauss,
          recenter;
  double  *buf,
          cos, sin,
          ax, ay,
          sigma_inc,
          sigma_mscale,
          *vec,
          *local_vec,
          fitrad,
          normrad,
          gain,
          x_orig, y_orig;
} PSF_STRUCT;

typedef struct
{
  int     nsec_x, nsec_y,
          nx, ny,
          hw,
          wdeg,
          ncomp,
          nvecs,
          ntot,
          nwxy,
          *deg;
  float   *sig;
  double  *chi2_n,
          **vec,
          **vecs;
} KER_STRUCT;

/*** END ***/
