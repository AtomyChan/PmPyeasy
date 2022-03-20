/* 
 * normalize - normalize an image to a mean value of 1.0
 *
 * USAGE: 
 *   normalize [-Bsc:ec,sr:er] infile outfile 
 *
 * OPTIONS:
 *   infile  the FITS file to be normalized
 *   outfile the normalized FITS file to create
 *   -Bsc:ec,sr:er = coordinates of box for normalizing, default: whole image
 *   -s = output short integer (BITPIX=16), default is FLOAT
 *   -l = output long integer (BITPIX=32), default is FLOAT
 *   -h = print "help" message (usage message) and exit
 *
 * Normalizes an image to a mean pixel value over the entire image of 1.0
 *
 * COMPILING:
 *
 * normalize uses the CCDPRoc.h header file and the same set of prototypes.
 * It is therefore able to be compiled using the gccfits script provided
 * with the PLANET CCDProc distribution:
 *
 *        gccfits normalize
 *
 * This guarantees that normalize is compiled with the HEASARC cfitsio
 * libraries.
 *
 * EXTERNAL LIBRARIES:
 *
 * normalize uses the HEASARC cfitsio library, which is available online 
 * at URL:
 *
 *        http://heasarc.gsfc.nasa.gov/docs/software/fitsio/
 *
 * AUTHOR:
 *
 *  R.W. Pogge, OSU Astronomy Dept.
 *  pogge@astronomy.ohio-state.edu
 *  1998 June 3
 *
 * Modification History:
 *   Version 0.0: 1998 Jun 3: based on imcombine.c, solves a problem
 *                            with IR flats that must be normalized
 *                            *after* on-off differencing.
 ***************************************************************************/

#include "CCDProc.h"

int 
main(int argc, char *argv[])
{
  
  fitsfile *infptr;     /* input FITS pointer 1 */ 
  fitsfile *outfptr;     /* output FITS pointer */

  char hist1[FLEN_VALUE];     /* history card */
  char datetime[FLEN_VALUE];  /* Formatted Date/Time Tag */
  int c;
  int status = 0;
  int nfound, anynull;
  float imgnorm = 1.0;        /* image normalization constant = 1/mean */
  long naxis1 = 2;
  long naxis2 = 2;
  long naxes1[2];
  long naxes2[2];
  long fpixel = 1;
  float nullval = 0.0; 
  long sc, ec, sr, er;        /* stats box coordinates */
  long npixels, i;
  float *indata;
  float *outdata;
  char infile[MED_STR_SIZE];        /* input filename */
  char outfile[MED_STR_SIZE];        /* output filename */

  int statbox = 0;          /* use whole image for stats by default */
  int bitpix = FLOAT_IMG;   /* output bitpix=-32 (float) by default */

  double bzero = 0.0;       /* default FITS image intensity scaling factors */
  double bscale = 1.0;  

  double mean, sigma;       /* mean and sigma of the image array (-n flag) */

  /* must be at least 2 arguments on the command line (+1 for command) */

  if (argc < 3) {
      Usage(); 
      exit(1); 
    }

  /* parse the command line, options first */

  while (--argc > 0 && (*++argv)[0]=='-')
    {
      c = *++argv[0];
      switch(c) {
      case 'h':             /* -h flag: print "help" and exit */
        Usage();
	exit(0);
        break;
      case 'B':             /* -B flag: coords of stats box */
	sscanf(*argv,"B%d:%d,%d:%d",&sc,&ec,&sr,&er);
	statbox = 1;
        break;
      case 'l':             /* -l flag: output BITPIX=32 */
        bitpix = LONG_IMG;
        break;
      case 's':             /* -s flag: output BITPIX=16 */
        bitpix = SHORT_IMG;
        break;
      default:
        fprintf(stderr,"normalize: illegal option %c\n",c);
        Usage();
        exit(1);
        break;
      }
    }

  /*
   * The first of the remaining arguments on the command line is the
   * name of the input image, the second is the output image.
   */

  strcpy(infile,argv[0]);
  strcpy(outfile,argv[1]);

  /* Open the input file and read it into memory */
        
  if ( fits_open_file(&infptr, infile, READONLY, &status) )
       FITSError(status,"Cannot open first FITS file");

  /* read the NAXIS1 and NAXIS2 keyword to get image size */

  if ( fits_read_keys_lng(infptr, "NAXIS", 1, 2, naxes1, &nfound, &status) )
       FITSError(status,"Cannot read first FITS file header");

  npixels  = naxes1[0] * naxes1[1];      /* number of pixels in the image */

  /* Allocate memory for the data as a 1-D floating arrays.  The first,
   * outdata, will be used for the first image in the stack, and as the
   * accumulation buffer for the summation.  The second will be the working
   * buffer for image 2-N in the list to be summed.
   */

  outdata = (float*) fvalloc(npixels);  
  indata = (float*) fvalloc(npixels);

  /* Read the first file in the list into outdata.  Keep it open since
   * we will eventually need its header as the template for the output
   * image. */

  if (fits_read_img(infptr, TFLOAT, fpixel, npixels, &nullval,
                    outdata, &anynull, &status))
    FITSError(status,"Cannot read input FITS file data");

  /* 
   * Normalize the image to 1.0 by computing the mean and sigma in
   * the stats box (whole array or -B flag arguments), and divide by
   * the derived normalization factor, imgnorm
   */

  if (statbox == 0) {  /* no statbox, use whole array */
    sc = 1;
    ec = naxes1[0];
    sr = 1;
    er = naxes1[1];
  }

  mean = immean(outdata,naxes1[0],naxes1[1],sc,ec,sr,er,&sigma);
    
  imgnorm = 1.0/(float)(mean);
  MulCon(outdata, imgnorm, npixels);

  /* make a history card with the mean & sigma used */

  TimeTag(datetime);
  sprintf(hist1,"(normalize) Mean=%.2f+/-%.2f",mean,sigma);

  /* compute the scaling parameters (BZERO & BSCALE) if doing integers */

  if (bitpix != FLOAT_IMG) {
    if (fitscl(outdata,npixels,bitpix,&bzero,&bscale))
      fprintf(stderr,"Error computing BZERO/BSCALE\n");
  }

  /* Write the results to a file or stdout */

  if (fits_create_file (&outfptr, outfile, &status))
    FITSError(status,"Cannot create output FITS file");

  /* copy the header from the input image */
  
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

  /* If doing integers, set the BZERO and BSCALE parameters */
  
  if (bitpix != FLOAT_IMG) {
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
    
  }
  
  /* write the data to the output image, datatype conversion is implicit */
  
  if (fits_close_file (infptr, &status))
    FITSError(status,"Cannot close input FITS file");

  if (fits_write_img (outfptr, TFLOAT, 1, npixels, outdata, &status))
    FITSError(status,"Cannot write output data"); 
   
  /* add the history card to the header */

  if (fits_write_history (outfptr, hist1, &status))
    FITSError(status,"Cannot update output history card(s)");

  /* Done! close the output FITS file */

  if (fits_close_file (outfptr, &status))
    FITSError(status,"Cannot close output FITS file");
  
  return 0;

}

void Usage()
{
  fprintf(stderr,"\nUsage: normalize [-Bslh] infile outfile\n");
  fprintf(stderr, " where infile  the FITS file to be normalized\n");
  fprintf(stderr, "       outfile the normalized FITS file to create\n");
  fprintf(stderr, "       -Bsc:ec,sr,er = coordinates of box for computing\n");
  fprintf(stderr, "            the normalization (-n flag)\n");
  fprintf(stderr, "            [Default: use whole image]\n");
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

void MulCon (float *img, float con, long n)
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

void TimeTag(char tagstr[FLEN_VALUE])
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

int fitscl(float *array, long npix, int bitpix, double *bzero, double *bscale)
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
 * immean - compute mean data value in a box
 *
 * Usage:
 *   immean(array,nc,nr,sc,ec,sr,er,niter,reject,>sigma)
 * Arguments:
 *   array (float): data array of size nc*nr (1-D array)
 *   nc, nr (long): dimensions of the input array
 *   sc,ec (long): starting & ending column pixels for the subarray (axis1)
 *   sr,er (long): starting & ending row pixels for the subarray (axis2)
 *   niter (long: number of interations [default: 1]
 *   reject (float): rejection threshold in units of sigma
 *   sigma (float): final computed sigma (standard deviation)
 *
 * If (sc,ec,sr,er) are 0, performs operation on the entire array.  If
 * niter=0, assumes 1 iteration (no sigma rejection)
 *
 * R. Pogge
 * OSU Astronomy
 * 1998 May 18
 *
 ***************************************************************************/

double immean(float *array, long nc, long nr, long sc, long ec, 
	      long sr, long er, double *sigma)
{

  long ic, ir, ii;
  long ipix, npix;
  double pixel;
  double mpix = 0.0;
  double sumx = 0.0;
  double sumxx = 0.0;
  long kount = 0;

  /* if the bounding box coords are 0, use the entire array */

  if (sc==0 && ec==0) {
    sc=1;
    ec=nc;
  }
  if (sr==0 && er==0) {
    sr=1;
    er=nr;
  }

  npix = (ec-sc+1)*(er-sr+1);

  if (npix < 2)
    return(0);

  /* Compute mean and sigma [the latter approximating for speed] */

  kount = 0;
  sumx = 0.0;
  sumxx = 0.0;
  for (ir = sr; ir < er+1; ir++)  {
    ii = (ir-1)*nc;
    for (ic = sc; ic < ec+1; ic++) {
      ipix = ii + ic - 1;
      pixel = (double) array[ipix];
      sumx = sumx + pixel;
      sumxx = sumxx + (pixel*pixel);
      kount++;
    }
  }
  mpix = (sumx)/(double)(npix);
  if (npix > 1) {
    *sigma = sqrt((sumxx-(sumx*sumx/(double)(npix)))/(double)(npix-1));
  } else {
    *sigma = 0.0;
  }

  return(mpix);
  
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

