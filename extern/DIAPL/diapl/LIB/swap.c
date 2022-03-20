/*========================================================*/
/*                                                        */
/*  swap.c           2003.03.04      version 0.1          */
/*                                                        */
/*  Copyright (C) 2003 Wojtek Pych, DDO                   */
/*                                                        */
/*  Written for GNU project C and C++ Compiler            */
/*                                                        */
/*========================================================*/

#include "swap.h"

/*--------------------------------------------------------*/
/* returns 1 on Intel architecture machines, 0 on Suns    */
/* Function taken from DIA package by P. Wozniak          */

int needswap(void)
{
        union { short s; char c[2]; } u;

  u.s=1;
  return((int)(u.c[0]));
}
/*--------------------------------------------------------*/
void swap_bytes(void *ptr, int ndata, int nbytes)
{
        char  c,
              *data;
        int   i, j, k, m;

  data=(char *)ptr;
  ndata*=nbytes;
  for (i=0; i<ndata; i+=nbytes)
  {
    k=nbytes/2;
    for (j=0; j<k; j++)
    {
      m=i+nbytes-1-j;
      c=data[i+j];
      data[i+j]=data[m];
      data[m]=c;
    }
  }

  return;
}

/*** END ***/
