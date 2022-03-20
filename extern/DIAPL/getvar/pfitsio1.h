/*========================================================*/
/*                                                        */
/*  pfitsio1.h      version 2.6.0         2006.10.25      */
/*                                                        */
/*  Copyright (C) 2006 by Wojtek Pych, CAMK PAN           */
/*                                                        */
/*  Written for GNU project C and C++ Compiler            */
/*  Private FITS I/O library.                             */
/*                                                        */
/*========================================================*/

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
float   *read_FITS_1Dfile(char *, char, int *, char ***, int *);
float   *read_FITS_2D1file(char *, char, int *, char ***, int *, int *);
void    scale(int, float *, short *);
void    write_FITS_1Dimage(FILE *, size_t, int, void *);
void    write_FITS_1Dfile(char *, int, char **, size_t, int, void *);
void    write_FITS_2D1file(char *, int, char **, size_t, int, void *);

/*** END ***/
