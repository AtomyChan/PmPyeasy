/*========================================================*/
/*                                                        */
/*  pfitsin.h       version 3.1.2         2005.11.08      */
/*                                                        */
/*  Copyright (C) 2005 by Wojtek Pych, CAMK PAN           */
/*  Copyright (C) 2004 by Wojtek Pych, DDO UofT           */
/*  Copyright (C) 2002 by Wojtek Pych, CAMK PAN           */
/*  Copyright (C) 1997-98 by Wojtek Pych, OAUW            */
/*                                                        */
/*  Written for GNU project C and C++ Compiler            */
/*  Private FITS I/O library.                             */
/*                                                        */
/*========================================================*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Some FITS magic numbers */
#define RECORD_SIZE   2880
#define RECORD_CARDS  36
#define CARD_SIZE     80
#define KEYWORD_SIZE  8
#define VALUE_SIZE    70

int   needswap(void);
void  swap_bytes(void *, int, int);
int   get_FITS_key(int, char **, char *, char*);
long  FITS_image_size(int, char **);
int   read_FITS_header(FILE *, char ***);
int   read_FITS_headers(FILE *, int **, char ****, long **);
int   read_FITS_image(FILE *, char **, int, void **);
float *read_FITS_1Dfile(char *, char, int *, char ***, int *);
float **read_FITS_2Dfile(char *, char, int *, char ***, int *, int *);

/*** END ***/
