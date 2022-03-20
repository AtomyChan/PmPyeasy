/*========================================================*/
/*                                                        */
/*  template.c          version 1.2.1   2005.04.29        */
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
/* Program makes a mosaic of fits images                  */
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
#include <string.h>

#include "errmess.h"
#include "pfitsio1.h"
#include "hedit.h"

#define NAME_LEN 200

int main(int argc, char *argv[])
{
        char  **newheader,
              **imnames,
              line[NAME_LEN],
              *inpfname, *outfname;
        int   nx0, ny0, nx, ny, k, l, i, j, nim, x_nim, y_nim, count, hsize;
        long  ofs, idx, first;
        double *im, *imout;
        FILE  *inpf;

/***/
  nx0 = ny0 = 0;
  imout = NULL;

/* IO stuff */

  if (argc != 5)
  {
    printf("\n\tUSAGE: template  image_names_file output_image x_nim y_nim\n");
    exit(-2);
  }

  inpfname = argv[1];
  outfname = argv[2];

  x_nim = atoi(argv[3]);
  y_nim = atoi(argv[4]);

  if (!(imnames=(char **)calloc(1, sizeof(char *))))
    errmess("calloc(imnames)");

  if (!(inpf=fopen(inpfname, "r"))) errmess(inpfname);

  for (i=0; !feof(inpf); i++)
  {
    fgets(line, NAME_LEN, inpf);
    if (feof(inpf)) break;

    if (!(imnames=(char **)realloc(imnames, (i+1)*sizeof(char *))))
      errmess("realloc(imnames)");

    if (!(imnames[i]=(char *)calloc(strlen(line)+1, sizeof(char))))
      errmess("calloc(imnames[i])");

    sscanf(line, "%s", imnames[i]);
  }

  fclose(inpf);
  nim=i;

  if (x_nim*y_nim != nim) errmess("wrong number of images");

  first = 1;
  count = 0;

/* loop over the images */

  for (l=0; l<y_nim; l++)
  {
    for (k=0; k<x_nim; k++)
    {
/* get test image */
      im=read_FITS_2D1file(imnames[count++], 's', &hsize, &newheader, &nx, &ny);

      if (first)
      {
        nx0 = nx*x_nim;
        ny0 = ny*y_nim;
        if (!(imout=(double *)malloc(nx0*ny0*sizeof(double))))
          errmess("malloc(imout)");
        first = 0;
      }

      idx = 0;
      ofs = nx0*ny*l + nx*k;

      for (j=0; j<ny; j++)
      {
        for (i=0; i<nx; i++)
        {
          imout[ofs++] = im[idx++];
        }
        ofs += nx0 - nx;
      }

      free(im);
    }
  }

  for (i=0; i<nim; i++) free(imnames[i]);
  free(imnames);

  if (count != nim) errmess("aftermath for number of images bad");

  hedit(hsize, newheader, "NAXIS1  ", 'i', "Axis length",
                         "                 %3d", nx0);
  hedit(hsize, newheader, "NAXIS2  ", 'i', "Axis length",
                         "                 %3d", ny0);

/* write output image */
  write_FITS_2D1file(outfname, hsize, newheader, nx0*ny0, sizeof(double),
                     (void *)imout);

  return(0);
}
