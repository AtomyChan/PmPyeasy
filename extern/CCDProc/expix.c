// 
// expix - extract pixel data with minimal header info
//
// Usage: expix [-s] infile outfile
//
// Keywords:
//    infile = input FITS file
//    outfile = output FITS file
//
//---------------------------------------------------------------------------

#include "CCDProc.h"

#define VERSION "expix v0.1 (2007-09-23)"

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

  // check for flags

  if (argc < 3) {
   Usage(); 
    exit(1); 
  }

  strcpy(infile,argv[1]);
  strcpy(outfile,argv[2]);

  printf("copying pixels from %s -> %s\n",infile,outfile);

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

  if (fits_close_file (infptr, &status))
    FITSError(status,"Cannot close input FITS file");

  // Open a FITS file for the output 

  if (fits_create_file (&outfptr, outfile, &status))
    FITSError(status,"Cannot create output FITS file");

  status = 0;
  fits_create_img(outfptr, bitpix, inaxis, inaxes, &status);
  if (status !=0) 
    FITSError(status,"Cannot create output FITS image structure"); 

  status = 0;
  fits_write_img (outfptr, TFLOAT, 1, npixels, indata, &status);
  if (status !=0) 
    FITSError(status,"Cannot write output FITS image"); 

  // Insert keyword ORIGFILE for the original filename

  status = 0;
  fits_update_key(outfptr,TSTRING,"ORIGFILE",infile,\
		  "Original Source File",&status);

  // Done! close the remaining open FITS files 

  status = 0;

  if (fits_close_file (outfptr, &status))
    FITSError(status,"Cannot close output FITS file");
  
  if (verbose) printf("expix done\n");

  return 0;

}

//*************************************************************************

void Usage()
{
  printf("\nUsage: expix inFITS outFITS\n");
  printf("  where inFITS = input FITS filename\n");
  printf("       outFITS = output FITS filename\n");
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

