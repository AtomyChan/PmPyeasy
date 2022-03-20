/*========================================================*/
/*                                                        */
/*  pfitsio.h       version 2.7.7         2006.10.19      */
/*                                                        */
/*  Copyright (C) 2006 by Wojtek Pych, CAMK PAN           */
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

int     needswap(void);
void    swap_bytes(void *, size_t, int);
int     get_FITS_key(int, char **, char *, char*);
size_t  FITS_image_size(int, char **);
int     read_FITS_header(FILE *, char ***);
int     read_FITS_headers(FILE *, int **, char ****, size_t **);
int     del_header_card(int, char ***, char *);
void    write_FITS_header(FILE *, int, char **);
int     rimage(FILE *, char **, int, void **);
size_t  read_FITS_image(FILE *, char **, int, void **);
double   *read_FITS_1Dfile(char *, char, int *, char ***, int *);
double   **read_FITS_2Dfile(char *, char, int *, char ***, int *, int *);
void    scale(int, double *, short *);
void    write_FITS_2Dimage(FILE *, int, int, int, void **);
void    write_FITS_2Dfile(char *, int, char **, int, int, int, void **);

/*** END ***/
