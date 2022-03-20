/*========================================================*/
/*                                                        */
/*  errmess.c       2006.10.19      version 2.1           */
/*                                                        */
/*  Copyright (C) 2006 by Wojtek Pych, CAMK PAN           */
/*                                                        */
/*  Written for GNU project C and C++ Compiler            */
/*                                                        */
/*  Print error message to stderr and terminate program.  */
/*                                                        */
/*========================================================*/

#include "errmess.h"

void errmess(char *errmes)
{
  fprintf(stderr, "\n\t\a");
  perror(errmes);
  fprintf(stderr, "\n");
  exit(EXIT_FAILURE);
}

/*** END ***/
