/*========================================================*/
/*                                                        */
/*  write_kernel.c      version 1.4.2   2005.04.20        */
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

#include "defs.h"
#include "errmess.h"

void write_kernel(char *fname, double *exp_vec, int nsec_x, int nsec_y,
                  int mode, PARAMS par, double *chi2_n)
{
        int   i;
        FILE  *outf;

  if (par.verbose > 2)
  {
    printf("Writing kernel to '%s'\n", fname);
    printf("mode= %d\n", mode);
  }

  if (mode == 0)
  {
    if (par.verbose > 4)
    {
      printf("nsec_x= %d   nsec_y= %d\n", nsec_x, nsec_y);
      printf("par.nx= %d   par.ny= %d\n", par.nx, par.ny);
      printf("par.kerhw= %d\n", par.kerhw);
      printf("par.wdeg= %d\n", par.wdeg);
      printf("par.ncomp= %d\n", par.ncomp);
      for (i=0; i<par.ncomp; i++)
        printf("  par.sig[%d]= %g\n", i, par.sig[i]);
      for (i=0; i<par.ncomp; i++)
        printf("  par.deg[%d]= %d\n", i, par.deg[i]);
      printf("--------\n");
    }

    if (!(outf=fopen(fname, "wb"))) errmess(fname);
    fwrite(&nsec_x    , sizeof(int), 1, outf);
    fwrite(&nsec_y    , sizeof(int), 1, outf);
    fwrite(&par.nx   , sizeof(int), 1, outf);
    fwrite(&par.ny   , sizeof(int), 1, outf);
    fwrite(&par.kerhw, sizeof(int), 1, outf);
    fwrite(&par.wdeg , sizeof(int), 1, outf);
    fwrite(&par.ncomp, sizeof(int), 1, outf);

    fwrite(par.sig, sizeof(double), par.ncomp, outf);
    fwrite(par.deg, sizeof(int)  , par.ncomp, outf);
  }
  else
  {
    if (!(outf=fopen(fname, "ab+"))) errmess(fname);
  }

  fwrite(chi2_n, sizeof(double), 1, outf);
  fwrite(&exp_vec[par.nbkg], sizeof(double), par.ntot-par.nbkg, outf);

  fclose(outf);

  if (par.verbose > 4)
  {
    printf("chi2_n= %g\n", *chi2_n);
    for (i=par.nbkg; i<par.ntot-par.nbkg; i++)
      printf("  exp_vec[%d]= %g\n", i, exp_vec[i]);
    printf("--------\n");
  }

  return;
}
