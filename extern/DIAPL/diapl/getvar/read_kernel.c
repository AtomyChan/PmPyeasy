/*========================================================*/
/*                                                        */
/*  read_kernel.c       version 1.3.0   2005.04.20        */
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
#include <stdlib.h>

#include "defs.h"
#include "errmess.h"

void read_kernel(char *kerfname, KER_STRUCT *ker, int mode, int verbose)
{
        int   nsec_x, nsec_y, wdeg, ncomp, ntot, nvecs, k, nwxy,
              nx, ny, kerhw, flag, *deg,
              i;
        double *sig;
        FILE  *inpf;

  if (verbose > 1) printf("read_kernel(%s)\n", kerfname);
  if (!(inpf = fopen(kerfname, "rb"))) errmess(kerfname);

  fread(&nsec_x , sizeof(int), 1, inpf);
  fread(&nsec_y , sizeof(int), 1, inpf);
  fread(&nx     , sizeof(int), 1, inpf);
  fread(&ny     , sizeof(int), 1, inpf);
  fread(&kerhw  , sizeof(int), 1, inpf);
  fread(&wdeg   , sizeof(int), 1, inpf);
  fread(&ncomp  , sizeof(int), 1, inpf);

  if (verbose > 2)
  {
    printf("nsec_x= %d\n", nsec_x);
    printf("nsec_y= %d\n", nsec_y);
    printf("nx    = %d\n", nx);
    printf("ny    = %d\n", ny);
    printf("kerhw = %d\n", kerhw);
    printf("wdeg  = %d\n", wdeg);
    printf("ncomp = %d\n", ncomp);
    printf("----------\n");
  }

  if (mode)
  {
    if (!(ker->sig = (double *)malloc(ncomp*sizeof(double))))
      errmess("malloc(ker->sig)");
    if (!(ker->deg = (int *)malloc(ncomp*sizeof(int))))
      errmess("malloc(ker->deg)");

    fread(ker->sig, sizeof(double), ncomp, inpf);
    fread(ker->deg, sizeof(int)  , ncomp, inpf);

    if (verbose > 2)
    {
      printf("mode= %d\n", mode);
      for (i=0; i<ncomp; i++)
        printf("  ker->sig[%d]= %g\n", i, ker->sig[i]);
      for (i=0; i<ncomp; i++)
        printf("  ker->deg[%d]= %d\n", i, ker->deg[i]);
      printf("----------\n");
    }

    nwxy = (wdeg+1)*(wdeg+2)/2;

    nvecs = 0;
    for (k=0; k<ncomp; k++)
      nvecs += (ker->deg[k]+1)*(ker->deg[k]+2)/2;

    ntot = nvecs + (nvecs - 1)*(nwxy - 1);

    if (verbose > 2) printf("nvecs= %d   ntot= %d\n", nvecs, ntot);

    if (!(ker->chi2_n = (double *)malloc(nsec_x*nsec_y*sizeof(double))))
      errmess("malloc(ker->chi2_n)");

    if (!(ker->vec    = (double **)malloc(nsec_x*nsec_y*sizeof(double *))))
      errmess("malloc(ker->vec)");
    for (k=0; k<nsec_x*nsec_y; k++)
      if (!(ker->vec[k] = (double *)malloc(ntot*sizeof(double))))
        errmess("malloc(ker->vec[k])");

    ker->nx    = nx;
    ker->ny    = ny;
    ker->hw    = kerhw;
    ker->wdeg  = wdeg;
    ker->nwxy  = nwxy;
    ker->ntot  = ntot;
    ker->ncomp = ncomp;
    ker->nvecs = nvecs;
    ker->nsec_x= nsec_x;
    ker->nsec_y= nsec_y;
  }
  else
  {
    flag = 0;

    if (ker->nsec_x != nsec_x)  flag=1;
    if (ker->nsec_y != nsec_y)  flag=1;
    if (ker->nx     != nx    )  flag=1;
    if (ker->ny     != ny    )  flag=1;
    if (ker->hw     != kerhw )  flag=1;
    if (ker->wdeg   != wdeg  )  flag=1;
    if (ker->ncomp  != ncomp )  flag=1;

    if (!(sig=(double *)malloc(ncomp*sizeof(double))))  errmess("malloc(sig)");
    if (!(deg=  (int *)malloc(ncomp*sizeof(int))))    errmess("malloc(deg)");

    fread(sig, sizeof(double), ncomp, inpf);
    fread(deg, sizeof(int)  , ncomp, inpf);

    for (k=0; k<ncomp; k++)
    {
      if (ker->deg[k] != deg[k])  flag = 1;
      if (ker->sig[k] != sig[k])  flag = 1;
    }

    free(deg);
    free(sig);

    if (flag)
    {
      printf("ERROR! read_kernel(): %s: parameters mismatch\n", kerfname);
      exit(5);
    }
  }

  for (k=0; k<nsec_x*nsec_y; k++)
  {
    fread(&ker->chi2_n[k], sizeof(double), 1, inpf);
    fread(ker->vec[k], sizeof(double), ker->ntot, inpf);
  }

  fclose(inpf);

  return;
}
