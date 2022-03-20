/* 
 * flt2int - convert BITPIX=-32 to BITPIX=16
 *
 * USAGE: 
 *   flt2int in32 out16
 *
 * OPTIONS:
 *   in32  the BITPIX=-32 (4-byte float) FITS file to be converted
 *   out16 the BITPIX=16 (2-byte integer) FITS file to create
 *   -h = print "help" message (usage message) and exit
 *
 * COMPILING:
 *
 * flt2int uses the CCDProc.h header file and the same set of prototypes.
 * It is therefore able to be compiled using the gccfits script provided
 * with the PLANET CCDProc distribution:
 *
 *        gccfits flt2int
 *
 * This guarantees that flt2int is compiled with the HEASARC cfitsio
 * libraries.
 *
 * EXTERNAL LIBRARIES:
 *
 * flt2int uses the HEASARC cfitsio library, which is available online 
 * at URL:
 *
 *        http://heasarc.gsfc.nasa.gov/docs/software/fitsio/
 *
 * AUTHOR:
 *
 *  R.W. Pogge, OSU Astronomy Dept.
 *  pogge@astronomy.ohio-state.edu
 *  2003 August 4
 *
 * Modification History:
 *
 ***************************************************************************/

#include "CCDProc.h"

int cleanbad(float *, long );

int 
main(int argc, char *argv[])
{
  
  fitsfile *infptr;     /* input FITS pointer 1 */ 
  fitsfile *outfptr;     /* output FITS pointer */

  int c;
  int status = 0;
  int nfound, anynull;
  long naxis1 = 2;
  long naxis2 = 2;
  long naxes1[2];
  long naxes2[2];
  long fpixel = 1;
  float nullval = 0.0; 
  long sc, ec, sr, er;        /* stats box coordinates */
  long npixels, i;
  float *imgdata;

  char infile[MED_STR_SIZE];        /* input filename */
  char outfile[MED_STR_SIZE];        /* output filename */

  int bitpix = SHORT_IMG;   /* output bitpix=16 (2-byte signed int) */

  double bzero = 0.0;       /* default FITS image intensity scaling factors */
  double bscale = 1.0;  

  /* must be at least 2 arguments on the command line (+1 for command) */

  if (argc < 3) {
      Usage(); 
      exit(1); 
    }

  /*
   * The first argument is the name of the input image, the second is
   * the output image.
   */

  strcpy(infile,argv[1]);
  strcpy(outfile,argv[2]);

  /* Open the input file and read it into memory */
        
  if ( fits_open_file(&infptr, infile, READONLY, &status) )
       FITSError(status,"Cannot open first FITS file");

  /* read the NAXIS1 and NAXIS2 keyword to get image size */

  if ( fits_read_keys_lng(infptr, "NAXIS", 1, 2, naxes1, &nfound, &status) )
       FITSError(status,"Cannot read first FITS file header");

  npixels  = naxes1[0] * naxes1[1];      /* number of pixels in the image */

  /*
   * Allocate memory for the data as a 1-D floating array.
   */

  imgdata = (float*) fvalloc(npixels);  

  /* 
   * Read the data into imgdata, but keep it open as we need the header later
   */

  if (fits_read_img(infptr, TFLOAT, fpixel, npixels, &nullval,
                    imgdata, &anynull, &status))
    FITSError(status,"Cannot read input FITS file data");

  /* DANGEROUS: clean out bad data */

  cleanbad(imgdata,npixels);

  /* compute the scaling parameters (BZERO & BSCALE) */

  if (fitscl(imgdata,npixels,bitpix,&bzero,&bscale))
    fprintf(stderr,"Error computing BZERO/BSCALE\n");

  /* Write the results to a file or stdout */

  if (fits_create_file (&outfptr, outfile, &status))
    FITSError(status,"Cannot create output FITS file");

  /* Copy the header from the input image */
  
  if (fits_copy_hdu (infptr, outfptr, 1, &status))
    FITSError(status,"Cannot copy header from first FITS file");

  /* Remove any delinquent BZERO or BSCALE cards (ignore errors) */

  fits_delete_key (outfptr,"BZERO",&status);
  status = 0;
  fits_delete_key (outfptr,"BSCALE",&status);
  status = 0;
  
  /* Modify the data type of the current array (use values for img 1) */
  
  if (fits_resize_img(outfptr,bitpix,naxis1,naxes1,&status))
    FITSError(status,"Cannot convert output data");

  /* Set the BZERO and BSCALE parameters */
  
  if (fits_set_bscale(outfptr,bscale, bzero, &status))
    FITSError(status,"Cannot scale output data"); 
    
  /* write the BZERO and BSCALE cards into the header (why CFITSIO doesn't
     do this when it scales we don't know) */

  if (fits_write_key_dbl(outfptr,"BSCALE",bscale,12,
			 "DATA=FITS*BSCALE+BZERO",&status))
    FITSError(status,"Cannot write BSCALE card"); 
    
  if (fits_write_key_dbl(outfptr,"BZERO",bzero,12,
			 "see BSCALE",&status))
    FITSError(status,"Cannot write BZERO card"); 
  
  /* write the data to the output image, datatype conversion is implicit */
  
  if (fits_write_img (outfptr, TFLOAT, 1, npixels, imgdata, &status))
    FITSError(status,"Cannot write output data"); 
   
  /* Done! close the output FITS file */

  if (fits_close_file (infptr, &status))
    FITSError(status,"Cannot close input FITS file");

  if (fits_close_file (outfptr, &status))
    FITSError(status,"Cannot close output FITS file");
  
  return 0;

}

void Usage()
{
  fprintf(stderr,"\nUsage: flt2int in32 out16\n");
  fprintf(stderr, " where in32  the BITPIX=-32 FITS file to convert\n");
  fprintf(stderr, "       out16 the BITPIX=16 FITS file to create\n");
}

/***************************************************************************
 * 
 * fvalloc - allocate virtual memory for a 1-D array of floats of size N
 *
 ***************************************************************************/

float* fvalloc(int N)
{
  float *out;
  
  if ((out = (float*) calloc(N,sizeof(float))) == NULL) {
    MEM_ERROR;
    exit(1);
  }
  
  return out;
}

/***************************************************************************
 *
 * cleanbad - clean out bad pixels
 *
 * Dangerous as hell, read carefully:
 *
 *    our presumption here is that the valid data in a single image run
 *    from 0 to 65534.  If this is not true, use another program.
 *
 *    All pixels with data values <0 are set to 0
 *    All pixels with data values >65534 are set to 65534
 *
 * You have been warned...
 *
 */

int 
cleanbad(float *array, long npix)
{

  long ii;
  float dmin, dmax;

  /* compute the data min and max values */

  dmin  = 0.0;
  dmax  = 65534.0;

  for (ii = 0; ii < npix; ii++)  {
    if ( array[ii] <= dmin )
      array[ii] = dmin;
    if ( array[ii] >= dmax )
      array[ii] = dmax;
  }

  return(0);
  
}

/***************************************************************************
 *
 * fitscl - compute FITS scaling parameters BZERO and BSCALE
 *
 * Arguments:
 *     array[npix] (i): floating array with data
 *     npix (i): number of pixels in the array
 *     bitpix (i): FITS bits/pixel (16 or 32)
 *     bzero (o): BZERO value
 *     bscale (o): BSCALE value
 *
 * 16 and 32-bit precision min/max are hardwired.  = 2**(BITPIX-1) - 1
 *
 * R. Pogge & P. Martini
 * OSU Astronomy
 * 1998 April 28
 *
 ***************************************************************************/

int 
fitscl(float *array, long npix, int bitpix, double *bzero, double *bscale)
{

  long ii;
  float dmin, dmax;
  double fitsmin, fitsmax;

  /* get data limits for bitpix.  If bitpix is not 16 or 32, error */

  if (bitpix == SHORT_IMG) {
    fitsmin = -32767.;
    fitsmax =  32767.;
  }
  else if (bitpix == LONG_IMG) {
    fitsmin = -2147483647.;
    fitsmax =  2147483647.;
  }
  else {
    fprintf(stderr,"Invalid BITPIX=%d, must be 16 or 32\n",bitpix);
    return(-1);
  }

  /* compute the data min and max values */

  dmin  = 1.0E30;
  dmax  = -1.0E30;

  for (ii = 0; ii < npix; ii++)  {
    if ( array[ii] < dmin )
      dmin = array[ii];
    
    if ( array[ii] > dmax )
      dmax = array[ii];
  }

  /* compute BZERO and BSCALE appropriate for BITPIX */

  *bzero = ((double)(dmin)*fitsmax-(double)(dmax)*fitsmin)/(fitsmax-fitsmin);
  if (dmax == dmin) {
    *bscale = 1.0;
  } else {
    *bscale = (double)(dmax-dmin) / (fitsmax-fitsmin);
  }
  return(0);
  
}

/*************************************************************************** 
 *
 * FITSError - Print out cfitsio error messages and exit program 
 *
 * Arguments:
 *     status (input, int): error status flag returned by the cfitsio
 *                          function
 *     msgstr (input, char): any additional information to print with
 *                           the error message (e.g., what it was doing
 *                           when it bombed to aid debugging)
 *
 * Translates status codes into text error messages and prints them
 * to stderr.  Since FITS errors are by default fatal, it exits the
 * program, passing the value of status to the shell.
 *
 * R. Pogge
 * OSU Astronomy Dept
 *
 * Modification History:
 *    1999 June 1: added msgstr argument [rwp/osu]
 *
 */ 

void FITSError(int status, char *msgstr)
{

  char logstr[MED_STR_SIZE];
  char status_str[FLEN_STATUS], fitserrmsg[FLEN_ERRMSG];
  int ierr = 0;

  if (status)
    fprintf(stderr,"Error: FITS error occurred during program execution\n");

  fprintf(stderr,"%s\n",msgstr);

  fits_get_errstatus(status, status_str);   /* get the error description */
  fprintf(stderr, "status=%d: %s\n", status, status_str);

  /* get first message; null if stack is empty */
  if ( fits_read_errmsg(fitserrmsg) )
    {
      fprintf(stderr, " %s\n", fitserrmsg);
      
      while ( fits_read_errmsg(fitserrmsg) )  { /* get remaining messages */
        fprintf(stderr, " %s\n", fitserrmsg);
      }
    }

  exit(status);       /* terminate the program, returning error status */

}

