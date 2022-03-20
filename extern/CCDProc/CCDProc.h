#ifndef CCDPROC_H
#define CCDPROC_H

/***************************************************************************
 *
 * CCDProc.h - header file for CCDProc package
 *
 * P. Martini & R. Pogge
 * OSU Astronomy Dept.
 * 1998 April 8
 *
 * Modification History:
 *   1998 May 20: added struct for systab & new error handling hooks [rwp/osu]
 *   1999 Aug 8: added dark-frame hooks to the systab struct [rwp/osu]
 *
 ***************************************************************************/

/******************************
 *
 * If SunOS, change the next line to read "#define __SUNOS"
 * if Solaris or Linux, change the next line to read "#undef __SUNOS"
 *
 * These are the only user-defined flags in this file 
 *
 ******************************/

#undef __SUNOS

/***************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/file.h>
#include <unistd.h>
#include <fcntl.h>   
#include <errno.h>
#include <time.h>
#include <sys/time.h>
#include "fitsio.h"  /* NASA/HEASARC CFITSIO Library header */

#define POSIX
#define BUF_SIZE        2048           /* generic buffer size */
#define LONG_STR_SIZE   4096           /* max size of a long string */
#define MED_STR_SIZE    256            /* max size of a medium string */
#define SHORT_STR_SIZE  32             /* max size of a short string */

#define ccd_TRUE  1
#define ccd_FALSE 0
#define SYSERR    -32768
#define NUL       '\0'

/* Some stock error handling stuff */

#define SIZE_MISMATCH fprintf(stderr, "\nError: Image sizes do not match\n")
#define MEM_ERROR     fprintf(stderr, "\nError: Out of memory\n")

#ifdef __SUNOS
int sys_nerr;
char *sys_errlist[];
#define ERRORSTR sys_errlist[errno]
#else
#define ERRORSTR strerror(errno)
#endif

/* Global system table structure   */

struct st {    
  int   fd_logfile;                  /* log file descriptor */
  int   dolog;                       /* keep a runtime processing log */
  int   bitpix;                      /* output FITS file data type */
  int   dobias;                      /* subtract overscan bias */
  int   namps;                       /* number of amps (max: 4)*/
  long  bsc[4],bec[4],bsr[4],ber[4]; /* bounding boxes of the overscan regions */
  long  dsc[4],dec[4],dsr[4],der[4]; /* bounding boxes of the data regions */
  int   dozero;                      /* subtract a zero (2D bias) frame */
  char  zerofile[MED_STR_SIZE];      /* zero (2D bias) image filename */
  int   dodark;                      /* subtract a dark frame */
  char  darkfile[MED_STR_SIZE];      /* dark frame filename */
  float darkscale;                   /* dark frame scaling factor [default=1] */
  int   doflat;                      /* divide by a flat-field frame */
  char  flatfile[MED_STR_SIZE];      /* flat-field image filename */
  int   dotrim;                      /* trim (crop) output image */
  long  tsc,tec,tsr,ter;             /* trim box coords */
  int   dofwhm;                      /* compute FWHM of stars & SKY level */
  float detgain;                     /* detector gain in e-/ADU */
  int   doclip;
  float minclip;
  float maxclip;
};

extern struct st *systab;   /* Global system table pointer   */

/* Function Prototypes */

void FITSError(int, char *);      /* Print FITSIO error messages */
void PrintError(char *);          /* Print Generic Error Messages */
void Usage();                     /* Print Usage statement for module */
float *fvalloc(int);              /* Allocate floating array of a given length */
void TimeTag(char *);                   /* Generate a date/time tag */
int fitscl(float *, long, int, double *, double *);  /* compute BZERO/BSCALE */

void ImAdd(float *, float *, long);  /* add two images together */
void ImSub(float *, float *, long);  /* subtract one image from another */
void ImMul(float *, float *, long);  /* multiply two images together */
void ImDiv(float *, float *, long);  /* divides one image by another */

void AddCon(float *, float , long);  /* add a constant to an image */
void AddBox (float *, int, int, int, int, int, float); /* add in box */

void ImClip(float *, float, float, long);

void MulCon(float *, float , long);  /* multiply an image by a constant */

int imtrim(float *,long,long,float *,long,long,long,long);  /* trim array */

double immean(float *,long,long,long,long,long,long,double *); /* image stats */

void ParseProc(char *);       /* process control file parser */

int  BCmp(char *, char *, int);     /* bitwise compare of 2 strings */
void GetArg(char *, int, char *);   /* get argument from line */
void UpperCase(char *);             /* convert string to upper case */
void BZero(char *, int);            /* clear (zero) a string */
void LeftStr(char *, char *, int);
void RightStr(char *, char *, int);
void MidStr(char *, char *, int, int);
void ltos(char *, long);

void gparam(int, int, float *, float, float *, float *, float *, float *);

#endif // CCDPROC_H
