/* 
 * imcombine - sum (or average) a set of FITS images
 *
 * USAGE: 
 *   imcombine [-anBslh] outfile infile1 infile2 ...
 *
 * OPTIONS:
 *   outfile is the FITS image to create
 *   infile1 infile2 ... are the input FITS images
 *   -a = average instead of sum
 *   -n = normalize to mean of 1.0
 *   -Bsc:ec,sr,er = coordinates of box for -n, default is whole image
 *   -s = output short integer (BITPIX=16), default is FLOAT
 *   -l = output long integer (BITPIX=32), default is FLOAT
 *   -h = do not copy header from the primary image [default: copy headers]
 *
 * Combines together a set of images by either summation or averaging
 * (simple unweighted arithmetic mean).  Can also generate an image
 * normalized to mean=1.0 by computing the mean data value within a given
 * region of interest.
 *
 * EXAMPLES:
 *
 *        imcombine -l mb9804.fits mb9804.001 mb9804.002 mb9804.003
 *
 * sum together mb9804.001 thru mb9804.003, creating the file mb9804.fits
 * with long-integer (BITPIX=32) format data.  Another way to do this would
 * have been to type:
 *
 *        imcombine -l mb9804.fits mb9804.00[1-3]
 *
 * or
 *
 *        imcombine -l mb9804.fits mb9804.00*
 *
 * if all of the mb9804.00* files in the current directory are the ones of
 * interest.  In general, use wildcards with caution.
 *
 * A mean of a stack of images may be computed using the command:
 *
 *        imcombine -a mzero.fits zero*.fits
 *
 * which creates a mean zero frame, mzero.fits, from zero frames zero*.fits
 * in the current working directory.
 *
 * To make a normalized mean flat, you would issue the following command:
 *
 *        imcombine -a -n -B17:528,1:512 mvflat.fits vflat.00*.fits
 *
 * which creates a mean normalized flat, mvflat.fits, from files
 * vflat.00*.fits, normalizing to 1.0 within the region of interest defined
 * by the box bounded by coordinates [17:528,1:512].
 *
 * COMPILING:
 *
 * imcombine uses the CCDPRoc.h header file and the same set of prototypes.
 * It is therefore able to be compiled using the gccfits script provided
 * with the PLANET CCDProc distribution:
 *
 *        gccfits imcombine
 *
 * This guarantees that imcombine is compiled with the HEASARC cfitsio
 * libraries.
 *
 * EXTERNAL LIBRARIES:
 *
 * imcombine uses the HEASARC cfitsio library, which is available online 
 * at URL:
 *
 *        http://heasarc.gsfc.nasa.gov/docs/software/fitsio/
 *
 * DISTRIBUTION:
 *
 * imcombine will be distributed with version v1.2 of the ccdproc
 * program, although it is currently available as an add-on with
 * v1.1 (current version at URL:
 *
 *        http://www.astronomy.ohio-state.edu/~pogge/PLANET/
 *
 * AUTHOR:
 *
 *        R.W. Pogge
 *        OSU Astronomy Dept.
 *        3 June 1998
 *
 * Modification History:
 *   Version 0.0: 1998 Jun 3: based on sumfits.c and mean.c, using the
 *                            current CCDProc prototypes. [rwp/osu]
 *
 *           0.1: 1998 Sep 2: cleaned up bugs, added -h flag [rwp/osu]
 *
 *   Version 1.0: 1998 Sep 28: released as an add-on to ccdproc [rwp/osu]
 *
 *           1.1: 2001 Oct 18: fixed bug in the -a feature [rwp/osu]
 *
 ***************************************************************************/

#include "CCDProc.h"

int 
main(int argc, char *argv[])
{
  
  fitsfile *in1fptr;     /* input FITS pointer 1 */ 
  fitsfile *in2fptr;     /* input FITS pointer 2 */
  fitsfile *outfptr;     /* output FITS pointer */

  char hist1[FLEN_VALUE];     /* history card */
  char hist2[FLEN_VALUE];     /* another history card */
  char datetime[FLEN_VALUE];  /* Formatted Date/Time Tag */
  int c;
  int status = 0;
  int nfound, anynull;
  int img;                    /* image counter */
  int nimgs;                  /* number of images to be summed together */
  float avgnorm = 1.0;        /* averaging normalization = 1/nimgs */
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
  char infile1[MED_STR_SIZE];        /* input file 1 */
  char infile2[MED_STR_SIZE];        /* input file 2-N */
  char outfile[MED_STR_SIZE];        /* output filename */

  int normalize = 0;        /* do not normalize the image to 1.0 by default */
  int average = 0;          /* sum rather than average images together */
  int statbox = 0;          /* use whole image for stats by default */
  int nohead = 0;           /* copy the primary image header by default */
  int bitpix = FLOAT_IMG;   /* output bitpix=-32 (float) by default */

  double bzero = 0.0;       /* default FITS image intensity scaling factors */
  double bscale = 1.0;  

  double mean, sigma;       /* mean and sigma of the image array (-n flag) */

  /* must be at least 3 arguments on the command line (+1 for command) */

  if (argc < 4) {
      Usage(); 
      exit(1); 
    }

  /* parse the command line, options first */

  while (--argc > 0 && (*++argv)[0]=='-')
    {
      c = *++argv[0];
      switch(c) {
      case 'a':             /* -a flag: average instead of sum */
        average = 1;
        break;
      case 'n':             /* -n flag: normalize image to 1.0 */
        normalize = 1;
        break;
      case 'l':             /* -l flag: output BITPIX=32 */
        bitpix = LONG_IMG;
        break;
      case 's':             /* -s flag: output BITPIX=16 */
        bitpix = SHORT_IMG;
        break;
      case 'h':             /* -h flag: print "help" and exit */
        nohead = 1;
        break;
      case 'B':             /* -B flag: coords of stats box */
	sscanf(*argv,"B%d:%d,%d:%d",&sc,&ec,&sr,&er);
	statbox = 1;
        break;
      default:
        fprintf(stderr,"imcombine: illegal option %c\n",c);
        Usage();
        exit(1);
        break;
      }
    }

  /* The first of the remaining arguments on the command line is the
   * name of the output image, the second is the first image of the list 
   * combined, up to the end of the command line.
   */

  nimgs = argc-1;
  if (nimgs < 2) {
    fprintf(stderr,"Must have at least 2 images to combine: nimgs=%d\n",nimgs);
    Usage();
    exit(1);
  }

  avgnorm = 1.0/(float)(nimgs);

  strcpy(outfile,argv[0]);
  strcpy(infile1,argv[1]);

  /* Open the input file and read it into memory */
        
  if ( fits_open_file(&in1fptr, infile1, READONLY, &status) )
       FITSError(status,"Cannot open first FITS file");

  /* read the NAXIS1 and NAXIS2 keyword to get image size */

  if ( fits_read_keys_lng(in1fptr, "NAXIS", 1, 2, naxes1, &nfound, &status) )
       FITSError(status,"Cannot read first FITS file header");

  npixels  = naxes1[0] * naxes1[1];      /* number of pixels in the image */

  /* Allocate memory for the data as a 1-D floating arrays.  The first,
   * outdata, will be used for the first image in the stack, and as the
   * accumulation buffer for the summation.  The second will be the working
   * buffer for image 2-N in the list to be summed.
   */

  outdata = (float*) fvalloc(npixels);  
  indata = (float*) fvalloc(npixels);

  /* Prepare a HISTORY card with the time/date and # of images */
  
  TimeTag(datetime);
  if (average==0) {
    sprintf(hist1,"%s:(imcombine) Sum of %d images",datetime,nimgs);
  } else {
    sprintf(hist1,"%s:(imcombine) Average of %d images",datetime,nimgs);
  }

  /* Read the first file in the list into outdata.  Keep it open since
   * we will eventually need its header as the template for the output
   * image. */

  if (fits_read_img(in1fptr, TFLOAT, fpixel, npixels, &nullval,
                    outdata, &anynull, &status))
    FITSError(status,"Cannot read first FITS file data");

  /* Read in each of the remaining images in turn, adding them to the
   * summation buffer */

  for (img=2; img<=nimgs; img++) {

    /* Open the next image on the list */

    if (fits_open_file(&in2fptr, argv[img], READONLY, &status))
      FITSError(status,"Cannot open current FITS image");

    /* read the NAXIS1 and NAXIS2 keyword to get image size.  If there
       is a size mismatch, abort now */
    
    if ( fits_read_keys_lng(in2fptr,"NAXIS",1,2,naxes2,&nfound,&status))
      FITSError(status,"Cannot read current FITS image header");

    if (naxes2[0]!=naxes1[0] || naxes2[1]!=naxes1[1]) {
      SIZE_MISMATCH; 
      fprintf(stderr,"Offending image: %s\n",argv[img]);
      exit(1); 
    }

    /* read the image into the working buffer */

    if (fits_read_img(in2fptr, TFLOAT, fpixel, npixels, &nullval,
		      indata, &anynull, &status))
      FITSError(status,"Cannot read current FITS image data");

    /* Add the working buffer to the summation buffer */

    ImAdd(outdata, indata, npixels);

    /* Close the current image and get the next on the list */

    if (fits_close_file (in2fptr, &status))
      FITSError(status,"Cannot close current FITS image data");

  }    

  /* If averaging, divide the summation buffer by the number of images
   * ("divide" = multiply by avgnorm=1/nimgs) 
   */

  if (average==1) 
    MulCon(outdata, avgnorm, npixels);

  /* 
   * If normalizing the image to 1.0, compute the mean and sigma in
   * the stats box (whole array or -B flag arguments), and divide by
   * the derived normalization factor, imgnorm
   */

  if (normalize==1) {

    if (statbox == 0) {  /* no statbox, use whole array */
      sc = 1;
      ec = naxes1[0];
      sr = 1;
      er = naxes1[1];
    }

    mean = immean(outdata,naxes1[0],naxes1[1],sc,ec,sr,er,&sigma);
    
    imgnorm = 1.0/(float)(mean);
    MulCon(outdata, imgnorm, npixels);

    /* make a second history card with the mean & sigma used */

    TimeTag(datetime);
    sprintf(hist2,"%s:(imcombine) Normalized with Mean=%.2f+/-%.2f",datetime,
	    mean,sigma);

  }

  /* compute the scaling parameters (BZERO & BSCALE) if doing integers */

  if (bitpix != FLOAT_IMG) {
    if (fitscl(outdata,npixels,bitpix,&bzero,&bscale))
      fprintf(stderr,"Error computing BZERO/BSCALE\n");
  }

  /* Write the results to a file or stdout */

  if (fits_create_file (&outfptr, outfile, &status))
    FITSError(status,"Cannot create output FITS file");

  /* copy the header from the input image */

  status = 0;

  if (nohead) {
    if (fits_create_img(outfptr, bitpix, naxis1, naxes1, &status))
      FITSError(status,"Cannot create primary header for output FITS file") ;
  } else {
    if (fits_copy_hdu (in1fptr, outfptr, 0, &status))
      FITSError(status,"Cannot copy header from first FITS file");
  }

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

  if (fits_close_file (in1fptr, &status))
    FITSError(status,"Cannot close first FITS file");

  if (fits_write_img (outfptr, TFLOAT, 1, npixels, outdata, &status))
    FITSError(status,"Cannot write output data"); 
   
  /* add the first history card */

  if (fits_write_history (outfptr, hist1, &status))
    FITSError(status,"Cannot update output history card(s)");

  /* add the second if normalizing */

  if (normalize==1) {
    if (fits_write_history (outfptr, hist2, &status))
      FITSError(status,"Cannot write normalization history");
  }

  /* Done! close the output FITS file */

  if (fits_close_file (outfptr, &status))
    FITSError(status,"Cannot close output FITS file");
  
  return 0;

}

void Usage()
{
  fprintf(stderr,"\nUsage: imcombine [-anBslh] outfile infile1 infile2 ...\n");
  fprintf(stderr, " where outfile is the FITS image to create\n");
  fprintf(stderr, "       infile1 infile2 ... are the input FITS images\n");
  fprintf(stderr, "       -a = compute the average instead of sum\n");
  fprintf(stderr, "            [Default: compute the sum]\n");
  fprintf(stderr, "       -n = normalize result to mean value of 1.0\n");
  fprintf(stderr, "            [Default: no normalization]\n");
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
 * ImAdd - add two images together
 *
 * arguments:
 *   im1 = image 1 (float)
 *   im2 = image 2 (float)
 *     n = total number of pixels in both images (long)
 *
 * returns:
 *   replaces im1 by im1+im2
 * 
 * notes:
 *   images must be perfectly aligned (same number of pixels and formats)
 *
 ***************************************************************************/

void ImAdd (float *im1, float *im2, long n)
{
  int i;
  register float *a, *b; 
  long N;
  
  a = im1; 
  b = im2; 
  N = n; 

  for (i=0; i<N; i++)
    *a++ += *b++; 

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

