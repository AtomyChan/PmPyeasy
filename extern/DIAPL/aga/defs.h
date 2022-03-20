/*========================================================*/
/*                                                        */
/*  defs.h              version 1.4.0   2005.04.19        */
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

#define MAX_NCOMP 5
#define MAX_FNAME 128
#define BIG_FLOAT 1.0e32
#define CHI2_INIT 999.0

#define VERB_NONE 0
#define VERB_LOW  1
#define VERB_MED  2
#define VERB_HIGH 3

#define RJCT_NONE 0
#define RJCT_LOW  1
#define RJCT_HIGH 2

typedef struct
{
  int     xc, yc,
          reject;
  double  chi2_n,
          **mat0,
          **mat1,
          **mat,
          *vec0,
          *vec;
} DOM_TABS;

typedef struct
{
  char  fits_ext[9];
  int   nx0, ny0,
        nx, ny,
        kerhw,
        ncomp,
        nvecs,
        bad_growrad1, bad_growrad2,
        bdeg,
        deg_inc,
        n_iter, 
        sdeg,
        wdeg,
        nbkg,
        nwxy,
        nspat,
        ntot,
        ndom,
        domhw,
        ndom_x, ndom_y,
        mohw,
        domain_mode,
        n_iter_dom,
        min_nkeep,
        *deg,
        verbose;
  float *sig,
        sig_inc,
        sat_level,
        min_level,
        gain,
        bad_value,
        min_area,
        n_sig,
        max_chi2,
        n_sig_dom,
        dom_thresh,
        min_area_dom;
} PARAMS;

/*** END ***/
