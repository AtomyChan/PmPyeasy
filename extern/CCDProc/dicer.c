// 
// dicer - dice an image into subimages
//
// Usage: dicer [-s] infile dicefile
//
// Keywords:
//    infile = input FITS file
//    dicefile = file with the dicing instructions
//
// If -s given, outputs data as short integers (BITPIX=16).  Default
// is BITPIX=-32 (floats)
//
//---------------------------------------------------------------------------

#include "CCDProc.h"

#define VERSION "dicer v0.1 (2004-09-22)"

int 
main(int argc, char *argv[])
{
  fitsfile *infptr;     // input image FITS file pointer  
  fitsfile *outfptr;    // output image FITS file pointer 
  
  int c;
  int status = 0;
  int nscan = 0;
  int nfound, anynull;
  long inaxis = 2;   // all images are 2-D FITS 
  long onaxis = 2;
  long inaxes[2];    // image dimensions 
  long onaxes[2];
  long fpixel = 1;
  long npixels, outpix, i, j;
  float *indata, *outdata;
  float nullval = 0.0; 
  int bitpix = FLOAT_IMG;              // default bitpix is -32 (float)

  long sc, ec, nc, sr, er, nr;

  char infile[MED_STR_SIZE] = "-"; // default input is STDIN 
  char procfile[MED_STR_SIZE];  // proc file 
  char outfile[MED_STR_SIZE];   // output filename

  double bzero = 0.0, bscale = 1.0;  // default scaling factors  
  int verbose = 0;

  char errmsg[MED_STR_SIZE];           /* generic error string */

  // procfile stuff

  char inbuf[BUF_SIZE];                // Generic input buffer 
  FILE *rfp;  // procfile pointer

  // check for flags

  if (argc < 3) {
   Usage(); 
    exit(1); 
  }

  while (--argc > 0 && (*++argv)[0]=='-') {
    c = *++argv[0];
    switch(c) {
    case 's':
    case 'S':
      bitpix = SHORT_IMG;
      break;
    case 'l':
    case 'L':
      bitpix = LONG_IMG;
      break;
    case 'v':
      verbose = 1;
      break;
    default:
      sprintf(errmsg,"dicer: illegal command-line option %c",c);
      PrintError(errmsg);
      Usage();
      exit(1);
      break;
    }
  } 

  strcpy(infile,argv[0]);
  strcpy(procfile,argv[1]);

  // Open the input FITS file READONLY                                

  if ( fits_open_file(&infptr, infile, READONLY, &status) )
       FITSError(status,"Cannot open input FITS file");

  // read the NAXIS1 and NAXIS2 keyword to get image size 

  if ( fits_read_keys_lng(infptr, "NAXIS", 1, 2, inaxes, &nfound, &status) )
       FITSError(status,"Cannot read input FITS header");

  npixels  = inaxes[0] * inaxes[1];      // number of pixels in the image 

  // Allocate memory for the data as a 1-D floating array. 

  indata = (float*) fvalloc(npixels);  

  // Read in the input file 

  if (fits_read_img(infptr, TFLOAT, fpixel, npixels, &nullval,
                    indata, &anynull, &status))
    FITSError(status,"Cannot read input FITS data");

  //**** Voodoo for SHORT, close then open file again [rwp/osu 1999Jun01 ****

  if (fits_close_file (infptr, &status))
    FITSError(status,"Cannot close input FITS file");

  if (fits_open_file(&infptr, infile, READONLY, &status) )
       FITSError(status,"Cannot open input FITS file");

  // Open the proc file and read through it, it contains our dicing
  // instructions

  if(!(rfp=fopen(procfile, "r"))) {
    fprintf(stderr,"Error: dicer unable to open %s\n",procfile);
    fprintf(stderr,"Reason: %s\n", ERRORSTR);
    exit(SYSERR);
  }

  // 
  // Loop through file reading 80 byte blocks at a time.  Skip 
  // Skip comments (#'s) and blank lines as you go.
  //
  //  Example format:
  //  # x_l   x_u     y_l     y_u     filename        dest_dir
  //  1       1074    1       1074    run11.blah.fits blah11/

  if (verbose)
    printf("Dicing up %s[1:%d,1:%d]...\n",infile,inaxes[0],inaxes[1]);

  while(fgets(inbuf, 80, rfp)) {
    if((inbuf[0]!='#') && (inbuf[0]!='\n') && inbuf[0]!=NUL) {
      inbuf[80] = NUL;

      nscan = sscanf(inbuf,"%d %d %d %d %s",&sc,&ec,&sr,&er,outfile);

      if (nscan != 5) {
	printf("Error: could not parse dicer file line '%s'\n",inbuf);
	printf("       incorrect format? (xs xe ys ye outfile) - dicer aborting\n");
	exit(1);
      }

      if (sc < 1 || sc > inaxes[0]) {
	printf("Error: requested xl=%d out of range, must be 1..%d\n",sc,inaxes[0]);
	printf("       offending dicer line '%s'\n",inbuf);
	exit (1);
      }

      if (ec < 1 || ec > inaxes[0]) {
	printf("Error: requested xu=%d out of range, must be 1..%d\n",ec,inaxes[0]);
	printf("       offending dicer line '%s'\n",inbuf);
	exit (1);
      }

      if (sr < 1 || sr > inaxes[1]) {
	printf("Error: requested yl=%d out of range, must be 1..%d\n",sr,inaxes[1]);
	printf("       offending dicer line '%s'\n",inbuf);
	exit (1);
      }

      if (er < 1 || er > inaxes[1]) {
	printf("Error: requested yu=%d out of range, must be 1..%d\n",er,inaxes[1]);
	printf("       offending dicer line '%s'\n",inbuf);
	exit (1);
      }

      nc = ec - sc + 1;
      nr = er - sr + 1;
      npixels = nc*nr;

      // Allocate memory for the output data

      outdata = (float*) fvalloc(npixels);  

      // extract the subframe from the bigger image

      if (imtrim(indata,inaxes[0],inaxes[1],outdata,sc,ec,sr,er)) {
	PrintError("Error: ccdproc could not trim image");
	exit(SYSERR);
      }

      if (bitpix != FLOAT_IMG) {
	if (fitscl(outdata,npixels,bitpix,&bzero,&bscale)) {
	  PrintError("Error: ccdproc could not compute BZERO/BSCALE");
	  exit(SYSERR);
	}
      }
      
      // Open a FITS file for the output 

      if (fits_create_file (&outfptr, outfile, &status))
	FITSError(status,"Cannot create output FITS file");

      // copy the header from the master image 

      if (fits_copy_hdu (infptr, outfptr, 32, &status))
	FITSError(status,"Cannot copy header to output FITS file");

      // Remove any delinquent BZERO or BSCALE cards (ignore errors) 
      
      fits_delete_key (outfptr,"BZERO",&status);
      status = 0;
      fits_delete_key (outfptr,"BSCALE",&status);
      status = 0;

      // Modify the data type of the current array 

      onaxes[0] = nc;
      onaxes[1] = nr;

      if (fits_resize_img(outfptr,bitpix,onaxis,onaxes,&status))
	FITSError(status,"Cannot convert data for output FITS file");

      // If doing integers, set the BZERO and BSCALE parameters 

      if (bitpix != FLOAT_IMG) {
	if (fits_set_bscale(outfptr,bscale, bzero, &status))
	  FITSError(status,"Cannot set output FITS header scaling params"); 
	
	// write the BZERO and BSCALE cards into the header (why CFITSIO doesn't
	// do this when it scales we don't know) 
	
	if (fits_write_key_dbl(outfptr,"BSCALE",bscale,14,
			       "DATA=FITS*BSCALE+BZERO",&status))
	  FITSError(status,"Cannot write BSCALE card"); 
	
	if (fits_write_key_dbl(outfptr,"BZERO",bzero,14,
			       "see BSCALE",&status))
	  FITSError(status,"Cannot write BZERO card"); 
    
      }

      // write the data to the output image, datatype conversion is implicit 
      // TFLOAT designates the *SOURCE* format, not the output format

      status = 0;

      fits_write_img (outfptr, TFLOAT, 1, npixels, outdata, &status);
      if (status !=0) 
	FITSError(status,"Cannot write trimmed output FITS image"); 

      // Done! close the remaining open FITS files 

      status = 0;

      if (verbose) 
	printf("   created subimage %s[%d:%d,%d:%d]\n",outfile,sc,ec,sr,er);

      if (fits_close_file (outfptr, &status))
	FITSError(status,"Cannot close output FITS file");
  
      // Free the vm for the output subimage

      free(outdata);

    }

  }

  // BOTTOM OF LOOP...

  if (fits_close_file (infptr, &status))
    FITSError(status,"Cannot close input FITS file");

  if (verbose) printf("dicer done\n");

  return 0;

}

//************************************************************************** 
//
//  PrintError - Print generic error messages w/o exiting.
//
  

void PrintError(char *errstr)
{
  char logstr[MED_STR_SIZE];
  int ierr = 0;

  fprintf(stderr,"%s\n",errstr);
}

//*************************************************************************

void Usage()
{
  printf("\nUsage: dicer [-s|l] [-v] inFITS dicefile\n");
  printf("  where inFITS  = input FITS filename\n");
  printf("        dicefile = the dicing instructions\n");
  printf("options: -s = output subimages as BITPIX=16 short ints\n");
  printf("         -l = output subimages as BITPIX=32 long ints\n");
  printf("         -v = verbose output while dicing\n");
}

//---------------------------------------------------------------------------

//
//  Allocate memory for a 1-dimensional array of floats
// 

float* fvalloc(int N)
{
  float *out;
  
  if ((out = (float*) calloc(N,sizeof(float))) == NULL) {
    MEM_ERROR;
    exit(1);
  }
  
  return out;
}

//---------------------------------------------------------------------------
//
// fitscl - compute FITS scaling parameters BZERO and BSCALE
//
// Arguments:
//     array[npix] (i): floating array with data
//     npix (i): number of pixels in the array
//     bitpix (i): FITS bits/pixel (16 or 32)
//     bzero (o): BZERO value
//     bscale (o): BSCALE value
//
// 16 and 32-bit precision min/max are hardwired to 1-2**(BITPIX-1) and
// 2**(BITPIX-1)-1, respectively.
//
// R. Pogge & P. Martini
// OSU Astronomy
// 1998 April 28
//

int 
fitscl(float *array, long npix, int bitpix, double *bzero, double *bscale)
{
  long ii;
  double dmin, dmax;
  double pixel;
  double fitsmin, fitsmax;

  // get data limits for bitpix.  If bitpix is not 16 or 32, error 

  if (bitpix == SHORT_IMG) {
    fitsmin = -32767;
    fitsmax =  32767;
  }
  else if (bitpix == LONG_IMG) {
    fitsmin = -2147483647;
    fitsmax =  2147483647;
  }
  else {
    fprintf(stderr,"(fitscl) Invalid BITPIX=%d, must be 16 or 32\n",bitpix);
    return(-1);
  }

  // compute the data min and max values 

  dmin  = 1.0E30;
  dmax  = -1.0E30;

  for (ii = 0; ii < npix; ii++)  {
    pixel = (double)(array[ii]);
    if ( pixel < dmin )
      dmin = pixel;
    
    if ( pixel > dmax )
      dmax = pixel;
  }

  // compute BZERO and BSCALE appropriate for BITPIX 

  *bzero = ((dmin*fitsmax)-(dmax*fitsmin))/(fitsmax-fitsmin);
  if (dmax == dmin) {
    *bscale = 1.0;
  } else {
    *bscale = (dmax-dmin) / (fitsmax-fitsmin);
  }
  return(0);
  
}

//---------------------------------------------------------------------------
//
// imtrim - trim an array down to a given size
//
// Arguments:
//   inarray (float): original data array of size nc*nr (1-D array)
//   nc, nr (long): dimensions of the input array
//   outarray (float): output subarray (1-D array of correct subarray size)
//   sc,ec (long): starting & ending column pixels for the subarray (axis1)
//   sr,er (long): starting & ending row pixels for the subarray (axis2)
//
// R. Pogge
// OSU Astronomy
// 1998 May 18
//

int 
imtrim(float *inarray, long nc, long nr, float *outarray, 
       long sc, long ec, long sr, long er)
{

  long ic, ir;
  long ipin, ipout;
  long nco, nro;
  long ii0, io0;

  nco = ec-sc+1;
  nro = er-sr+1;

  for (ir = sr; ir < er+1; ir++)  {
    ii0 = (ir-1)*nc;
    io0 = (ir-sr)*nco;
    for (ic = sc; ic < ec+1; ic++) {
      ipin = ii0 + ic - 1;
      ipout = io0 + (ic-sc+1) - 1;
      outarray[ipout] = inarray[ipin];
    }
  }

  return(0);
  
}

//---------------------------------------------------------------------------
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

  if ( fits_read_errmsg(fitserrmsg) ) {
    fprintf(stderr, " %s\n", fitserrmsg);
    while ( fits_read_errmsg(fitserrmsg) )  { // get remaining messages 
      fprintf(stderr, " %s\n", fitserrmsg);
    }
  }

  exit(status);       // terminate the program, returning error status 

}

