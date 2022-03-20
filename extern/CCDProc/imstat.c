// 
// imstat - compute statistics for an image
//
// USAGE: 
//   imstat [-Bsc:ec,sr:er] infile
//
// OPTIONS:
//   infile  the FITS file to be analyzed
//   -Bsc:ec,sr:er = coordinates of box for analysis, default: whole image
//   -h = print "help" message (usage message) and exit
//
// Computes statistics on an image or image subsection, and reports
// the following data to stdout
//    min
//    max
//    mean
//    sigma
//    mode
//
// imstat uses the CCDPRoc.h header file and the same set of prototypes.
// It is therefore able to be compiled using the gccfits script provided
// with the PLANET CCDProc distribution:
//
//        gccfits imstat
//
// This guarantees that imstat is compiled with the HEASARC cfitsio
// libraries.
//
// EXTERNAL LIBRARIES:
//
// imstat uses the HEASARC cfitsio library, which is available online 
// at URL:
//
//        http://heasarc.gsfc.nasa.gov/docs/software/fitsio/
//
// AUTHOR:
//
//  R.W. Pogge, OSU Astronomy Dept.
//  pogge@astronomy.ohio-state.edu
//  1998 June 3
//
// Modification History:
//   Version 0.0: 2002 Jun 3: based on normalize.c, uses the sky stuff
//                            from ccdproc [rwp/osu]
//
//*************************************************************************

#include "CCDProc.h"

// function prototype

void imstat(float *, long , long , long , long , long , long , 
	    double *, double *, double *, double *, double *);

// Main Program

int 
main(int argc, char *argv[])
{
  
  fitsfile *infptr;     // input FITS pointer  

  int c;
  int status = 0;
  int nfound, anynull;
  long naxis1 = 2;
  long naxis2 = 2;
  long naxes1[2];
  long naxes2[2];
  long fpixel = 1;
  float nullval = 0.0; 
  long sc, ec, sr, er;        // stats box coordinates 
  long npixels, i;
  float *indata;
  char infile[MED_STR_SIZE];  // input filename 

  int statbox = 0;          // use whole image for stats by default 
  int bitpix = FLOAT_IMG;   // output bitpix=-32 (float) by default 

  double bzero = 0.0;       // default FITS image intensity scaling factors 
  double bscale = 1.0;  

  double mean, sigma;   // mean and sigma of the image array (-n flag) 
  double dmin, dmax;    // min & max data values
  double mode;          // mode

  // There must be at least 1 arguments on the command line (+1 for command) 

  if (argc < 2) {
    Usage(); 
    exit(1); 
  }

  // Parse the command line, options first 

  while (--argc > 0 && (*++argv)[0]=='-')
    {
      c = *++argv[0];
      switch(c) {
      case 'h':             // -h flag: print "help" and exit 
        Usage();
	exit(0);
        break;
      case 'B':             // -B flag: coords of stats box 
	sscanf(*argv,"B%d:%d,%d:%d",&sc,&ec,&sr,&er);
	statbox = 1;
        break;
      default:
        fprintf(stderr,"imstat: illegal option %c\n",c);
        Usage();
        exit(1);
        break;
      }
    }

  // The remaining arguments on the command line should be the name of
  // the input image
 
  strcpy(infile,argv[0]);

  // Open the input file and read it into memory 
        
  if ( fits_open_file(&infptr, infile, READONLY, &status) )
       FITSError(status,"Cannot open first FITS file");

  // Read the NAXIS1 and NAXIS2 keyword to get image size 

  if ( fits_read_keys_lng(infptr, "NAXIS", 1, 2, naxes1, &nfound, &status) )
       FITSError(status,"Cannot read first FITS file header");

  npixels  = naxes1[0] * naxes1[1];      // number of pixels in the image 

  // Allocate memory for the data as a 1-D floating array to hold the image
   
  indata = (float*) fvalloc(npixels);  

  // Read the FITS file in the list into indata and close it.

  if (fits_read_img(infptr, TFLOAT, fpixel, npixels, &nullval,
                    indata, &anynull, &status))
    FITSError(status,"Cannot read input FITS file data");

  if (fits_close_file (infptr, &status))
    FITSError(status,"Cannot close input FITS file");

  // Are we using a stats box or the whole image?

  if (statbox == 0) {
    sc = 1;
    ec = naxes1[0];
    sr = 1;
    er = naxes1[1];
  }

  // Compute basic image statistics

  imstat(indata,naxes1[0],naxes1[1],sc,ec,sr,er,
	 &dmin,&dmax,&mean,&mode,&sigma);

  // For now report data to stdout

  printf("%.2f %.2f %.3f %.3f %.3f\n",dmin,dmax,mean,mode,sigma);

  return 0;

}

void Usage()
{
  fprintf(stderr,"\nUsage: imstat [-Bh] infile\n");
  fprintf(stderr," where infile  the FITS file to be analyzed\n");
  fprintf(stderr,"       -Bsc:ec,sr:er = coordinates of box for analysis\n");
  fprintf(stderr,"            the normalization (-n flag)\n");
  fprintf(stderr,"            [Default: use whole image]\n");
  fprintf(stderr,"       -h = print this message and exit\n\n");
  fprintf(stderr,"Output:\n");
  fprintf(stderr,"   min max mean mode sigma\n\n");
}

//**************************************************************************
// 
// fvalloc - allocate virtual memory for a 1-D array of floats of size N
//
//*************************************************************************

float* fvalloc(int N)
{
  float *out;
  
  if ((out = (float*) calloc(N,sizeof(float))) == NULL) {
    MEM_ERROR;
    exit(1);
  }
  
  return out;
}

//**************************************************************************
//
// imstat - compute min/max/mean/mode/sigma data value in a box
//
// Arguments:
//   array (float): data array of size nc*nr (1-D array)
//   nc, nr (long): dimensions of the input array
//   sc,ec (long): starting & ending column pixels for the subarray (axis1)
//   sr,er (long): starting & ending row pixels for the subarray (axis2)
//   niter (long: number of interations [default: 1]
//   reject (float): rejection threshold in units of sigma
//   dmin (double): minimum data value
//   dmax (double): maximum data value
//   mean (double): mean data value 
//   mode (double): mode of the data distribution
//   sigma (double): sigma (standard deviation) about the mean
//
// If (sc,ec,sr,er) are 0, performs operation on the entire array.  If
// niter=0, assumes 1 iteration (no sigma rejection)
//
// The mode calculation simply populates a histogram for data values 
// between -1000 and 63000 (bot and top, respectively below), and finds
// the peak value.  Not the most robust algorithm, but it'll get you
// close.
//
// R. Pogge, OSU Astronomy
// 2002 May 18
//
//*************************************************************************

void
imstat(float *array, long nc, long nr, long sc, long ec, long sr, long er, 
       double *dmin, double *dmax, double *mean, double *mode, double *sigma)
{

  long ic, ir, ii;
  long ipix, npix;
  double pixel;
  double mpix = 0.0;
  double sumx = 0.0;
  double sumxx = 0.0;
  long kount = 0;

  double bot = -1000.0;  // minimum data value to allow in a histogram
  double top = 64000.0;  // maximum data value to allow in a histogram
  int nhist = 65000;  // width of the data histogram in number of bins
  int *hist;
  double aa, bb;
  int mx, i, jx, kk;
  int z;

  // If the bounding box coords are 0, use the entire array 

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
    return;

  // initialize variables

  kount = 0;
  sumx = 0.0;
  sumxx = 0.0;
  *dmin = 1.0E30;
  *dmax = -1.0E30;
  *mean = 0.0;
  *sigma = 0.0;

  // histogram normalizations

  aa = (double)(nhist) / (top-bot);
  bb = -aa * bot;
  hist = (int *)calloc(nhist,sizeof(int));  // space for the intensity histogram

  // Compute statistics, using an approximate computation of sigma for speed

  for (ir = sr; ir < er+1; ir++)  {
    ii = (ir-1)*nc;
    for (ic = sc; ic < ec+1; ic++) {
      ipix = ii + ic - 1;
      pixel = (double) array[ipix];
      
      // Test the pixel against the min/max values

      if (pixel < *dmin)
	*dmin = pixel;
      if (pixel > *dmax)
	*dmax = pixel;

      // populate the histogram

      if ((pixel < top) && (pixel > bot)) {
	kk = (int)(aa * pixel + bb);
	*(hist+kk) += 1;
      }

      // Sum(x) and Sum(x*x) buffers

      sumx = sumx + pixel;
      sumxx = sumxx + (pixel*pixel);

      // increment the pixel counter

      kount++;
    }
  }
  mpix = (sumx)/(double)(npix);
  if (npix > 1) {
    *sigma = sqrt((sumxx-(sumx*sumx/(double)(npix)))/(double)(npix-1));
  } else {
    *sigma = 0.0;
  }
  *mean = mpix;

  // find the mode of the histogram

  mx = 0;
  for (i=0; i<nhist; i++)
    if (*(hist+i)>mx) {
      jx = i;
      mx = *(hist+i);
    }

  *mode = 1 + ((jx - bb) / aa);

  // free the memory allocated for the history array

  free(hist);
  
}

//************************************************************************** 
//
// FITSError - Print out cfitsio error messages and exit program 
//
// Arguments:
//     status (input, int): error status flag returned by the cfitsio
//                          function
//     msgstr (input, char): any additional information to print with
//                           the error message (e.g., what it was doing
//                           when it bombed to aid debugging)
//
// Translates status codes into text error messages and prints them
// to stderr.  Since FITS errors are by default fatal, it exits the
// program, passing the value of status to the shell.
//
// R. Pogge
// OSU Astronomy Dept
//
// Modification History:
//    1999 June 1: added msgstr argument [rwp/osu]
//
  
void 
FITSError(int status, char *msgstr)
{

  char logstr[MED_STR_SIZE];
  char status_str[FLEN_STATUS], fitserrmsg[FLEN_ERRMSG];
  int ierr = 0;

  if (status)
    fprintf(stderr,"Error: FITS error occurred during program execution\n");

  fprintf(stderr,"%s\n",msgstr);

  fits_get_errstatus(status, status_str);   // get the error description 
  fprintf(stderr, "status=%d: %s\n", status, status_str);

  // get first message; null if stack is empty 
  if ( fits_read_errmsg(fitserrmsg) )
    {
      fprintf(stderr, " %s\n", fitserrmsg);
      
      while ( fits_read_errmsg(fitserrmsg) )  { // get remaining messages 
        fprintf(stderr, " %s\n", fitserrmsg);
      }
    }

  exit(status);       // terminate the program, returning error status 

}

