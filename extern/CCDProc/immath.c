/* 
 * immath - add or multiply an image by a constant
 *
 * USAGE: 
 *   immath [-a###|-m###] infile outFile 
 *
 * OPTIONS:
 *   infile  input FITS file
 *   outFile output FITS file
 *   -a### = add constant value ###
 *   -m### = multiply by value ###
 *   -s = output short integer (BITPIX=16), default is FLOAT
 *   -l = output long integer (BITPIX=32), default is FLOAT
 *   -h = print "help" message (usage message) and exit
 *
 * Adds or multiplies all pixels in an image by a given constant value.
 *
 * COMPILING:
 *
 * immath uses the CCDPRoc.h header file and the same set of prototypes.
 * It is therefore able to be compiled using the gccfits script provided
 * with the MicroFUN CCDProc distribution:
 *
 *        gccfits immath
 *
 * This guarantees that immath is compiled with the HEASARC cfitsio
 * libraries.
 *
 * EXTERNAL LIBRARIES:
 *
 * immath uses the HEASARC cfitsio library, which is available online 
 * at URL:
 *
 *        http://heasarc.gsfc.nasa.gov/docs/software/fitsio/
 *
 * AUTHOR:
 *
 *  R.W. Pogge, OSU Astronomy Dept.
 *  pogge@astronomy.ohio-state.edu
 *  2009 June 19
 *
 * Modification History:
 *   Version 0.0: 2009 Jun 19: surprised we didn't need this sooner...
 *
 ***************************************************************************/

#include "CCDProc.h"

int 
main(int argc, char *argv[])
{
  
  fitsfile *inPtr;     // input FITS pointer
  fitsfile *outPtr;    // output FITS pointer

  char hist1[FLEN_VALUE];     /* history card */
  char dateTime[FLEN_VALUE];  /* Formatted Date/Time Tag */
  int c;
  int status = 0;
  int numFound, anyNull;

  int doAdd = 0;
  int doMult = 0;
  float constVal = 0.0;  // constant value to apply

  long naxis1 = 2;
  long naxis2 = 2;
  long naxes1[2];
  long naxes2[2];
  long fpixel = 1;
  float nullVal = 0.0; 
  long numPix, i;
  float *inData;
  float *outData;
  char inFile[MED_STR_SIZE];        /* input filename */
  char outFile[MED_STR_SIZE];        /* output filename */

  int bitpix = FLOAT_IMG;   /* output bitpix=-32 (float) by default */

  double bzero = 0.0;       /* default FITS image intensity scaling factors */
  double bscale = 1.0;  

  // Must be at least 3 arguments on the command line (+1 for command)

  if (argc < 4) {
    Usage(); 
    exit(1); 
  }

  // Parse the command line, options first

  while (--argc > 0 && (*++argv)[0]=='-') {
    c = *++argv[0];
    switch(c) {
    case 'h':             // -h flag: print "help" and exit
      Usage();
      exit(0);
      break;
    case 'a':             // -a: add constant
      sscanf(*argv,"a%f",&constVal);
      doAdd = 1;
      doMult = 0;
      break;
    case 'm':             // -m: multiply by constant
      sscanf(*argv,"m%f",&constVal);
      doMult = 1;
      doAdd = 0;
      break;
    case 'l':             // -l flag: output BITPIX=32
      bitpix = LONG_IMG;
      break;
    case 's':             // -s flag: output BITPIX=16
      bitpix = SHORT_IMG;
      break;
    default:
      fprintf(stderr,"immath: illegal option %c\n",c);
      Usage();
      exit(1);
      break;
    }
  }

  if (doAdd == 0 && doMult == 0) {
    printf("ERROR: no operation given: must specify -a### or -m###\n");
    Usage();
    exit(1);
  }

  // The first of the remaining arguments on the command line is the
  // name of the input image, the second is the output image.

  strcpy(inFile,argv[0]);
  strcpy(outFile,argv[1]);

  // Open the input file and read it into memory
        
  if ( fits_open_file(&inPtr, inFile, READONLY, &status) )
       FITSError(status,"Cannot open first FITS file");

  // Read the NAXIS1 and NAXIS2 keyword to get image size 

  if ( fits_read_keys_lng(inPtr, "NAXIS", 1, 2, naxes1, &numFound, &status) )
       FITSError(status,"Cannot read first FITS file header");

  numPix  = naxes1[0] * naxes1[1];      // number of pixels in the image

  // Allocate memory for the data as a 1-D floating arrays.  The
  // first, outData, will be used for the first image in the stack,
  // and as the accumulation buffer for the summation.  The second
  // will be the working buffer for image 2-N in the list to be
  // summed.
   
  outData = (float*) fvalloc(numPix);  
  inData = (float*) fvalloc(numPix);

  // Read inFile into outData.  Keep it open as we still need its HDU

  if (fits_read_img(inPtr, TFLOAT, fpixel, numPix, &nullVal,
                    outData, &anyNull, &status))
    FITSError(status,"Cannot read input FITS file data");

  if (doAdd)
    AddCon(outData,constVal,numPix);
  else if (doMult)
    MulCon(outData,constVal,numPix);

  // make a history card with the mean & sigma used

  TimeTag(dateTime);
  if (doAdd)
    sprintf(hist1,"(immath) Added %f",constVal);
  if (doMult)
    sprintf(hist1,"(immath) Multiplied by %f",constVal);

  // Compute the scaling parameters (BZERO & BSCALE) if doing integers

  if (bitpix != FLOAT_IMG) {
    if (fitscl(outData,numPix,bitpix,&bzero,&bscale))
      fprintf(stderr,"Error computing BZERO/BSCALE\n");
  }

  // Write the results to the file

  if (fits_create_file (&outPtr, outFile, &status))
    FITSError(status,"Cannot create output FITS file");

  // copy the header from the input image
  
  if (fits_copy_hdu (inPtr, outPtr, 1, &status))
    FITSError(status,"Cannot copy header from first FITS file");

  // Remove any delinquent BZERO or BSCALE cards (ignore errors)

  fits_delete_key (outPtr,"BZERO",&status);
  status = 0;
  fits_delete_key (outPtr,"BSCALE",&status);
  status = 0;
  
  // Modify the data type of the current array
  
  if (fits_resize_img(outPtr,bitpix,naxis1,naxes1,&status))
    FITSError(status,"Cannot convert output data");

  // If doing integers, set the BZERO and BSCALE parameters
  
  if (bitpix != FLOAT_IMG) {
    if (fits_set_bscale(outPtr,bscale, bzero, &status))
      FITSError(status,"Cannot scale output data"); 
    
    // Write the BZERO and BSCALE cards into the header (why CFITSIO
    // doesn't do this when it scales we don't know)

    if (fits_write_key_dbl(outPtr,"BSCALE",bscale,12,
			   "DATA=FITS*BSCALE+BZERO",&status))
      FITSError(status,"Cannot write BSCALE card"); 
    
    if (fits_write_key_dbl(outPtr,"BZERO",bzero,12,
			   "see BSCALE",&status))
      FITSError(status,"Cannot write BZERO card"); 
    
  }
  
  // Write the data to the output image, datatype conversion is implicit
  
  if (fits_close_file (inPtr, &status))
    FITSError(status,"Cannot close input FITS file");

  if (fits_write_img (outPtr, TFLOAT, 1, numPix, outData, &status))
    FITSError(status,"Cannot write output data"); 
   
  // Add the history card to the header

  if (fits_write_history (outPtr, hist1, &status))
    FITSError(status,"Cannot update output history card(s)");

  // Done! close the output FITS file

  if (fits_close_file (outPtr, &status))
    FITSError(status,"Cannot close output FITS file");
  
  return 0;

}

void 
Usage()
{
  fprintf(stderr,"\nUsage: immath [-Bslh] inFile outFile\n");
  fprintf(stderr, " where inFile  the input FITS file\n");
  fprintf(stderr, "       outFile the output FITS file to create\n");
  fprintf(stderr, "       -a### = add ### to all pixels\n");
  fprintf(stderr, "       -m### = multiply by ### to all pixels\n");
  fprintf(stderr, "       -s = output short integer (BITPIX=16)\n");
  fprintf(stderr, "       -l = output long integer (BITPIX=32)\n");
  fprintf(stderr, "            [Default: FLOAT, BITPIX=-32]\n");
  fprintf(stderr, "       -h = print this message and exit\n\n");
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
 * MulCon - Multiply an image by a constant
 *
 * arguments:
 *   img = image (float)
 *   con = constant (float)
 *     n = total number of pixels in the image (long)
 *
 * returns:
 *   replaces each pixel in img by img[i,j]*con
 * 
 ***************************************************************************/

void 
MulCon (float *img, float con, long n)
{
  int i;
  register float *a, c; 
  long N;
  
  a = img; 
  c = con; 
  N = n; 

  for (i=0; i<N; i++)
    *a++ *= c; 

}

/***************************************************************************
 *
 * AddCon - Add a constant to an image
 *
 * arguments:
 *   img = image (float)
 *   con = constant (float)
 *     n = total number of pixels in the image (long)
 *
 * returns:
 *   replaces each pixel in img by img[i,j]+con
 * 
 ***************************************************************************/

void 
AddCon (float *img, float con, long n)
{
  int i;
  register float *a, c; 
  long N;
  
  a = img; 
  c = con; 
  N = n; 

  for (i=0; i<N; i++)
    *a++ += c; 

}


/***************************************************************************
 *
 * TimeTag() - Create a date/time string by reading the system clock
 *
 * arguments
 *     tagstr - string to contain the timetag.
 *
 * Reads the local time & uses strftime() to format the output string.
 * FLEN_VALUE comes from fitsio.h
 *
 * R. Pogge
 * OSU Astronomy Dept.
 * 1998 April 9
 *
 ***************************************************************************/

void 
TimeTag(char tagstr[FLEN_VALUE])
{
  time_t now;
  struct tm *tm;

  now = time(NULL);
  tm = localtime(&now);
  if (strftime(tagstr,FLEN_VALUE,"%c",tm)==0)
    sprintf(tagstr,"Time Unavailable");
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

void 
FITSError(int status, char *msgstr)
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

