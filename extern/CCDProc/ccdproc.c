/* 
 * ccdproc - process CCD images
 *
 * Processes CCD images according to the parameters in an external
 * process control file.  Operations include, in this order:
 *       subtract overscan bias (up to 4 amps)
 *       trim the image (to exclude overscan)
 *       subtract a zero (2D bias) frame
 *       divide by a normalized flat field frame
 *       estimate DoPhot input parameters (sky, star fwhm, & thresholds)
 *
 * Usage: ccdproc -Pprocfile infile outfile
 *
 * Keywords:
 *    -Pprocfile = process control file.  See ccdproc.hlp for processing 
 *                 and a detailed discussion of the operation of this program
 *    infile = input FITS file
 *    outfile = output FITS file
 *
 * ccdproc does not permit overwriting of the input image.
 *
 * See ccdproc.html for complete documentation.
 *
 * Rick Pogge (pogge@astronomy.ohio-state.edu)
 * Paul Martini (martini@astronomy.ohio-state.edu)
 * OSU Astronomy Department
 * 1998 May 19
 *
 * Modification History:
 * 0.0 1998 Mar 19: began mods of PM's fccdproc [rwp/osu]
 * 0.1 1998 May 19: bulldozed and added procfile function [rwp/osu]
 * 0.2 1998 May 22: added fwhmsky function following code by John Menzies 
 *                  at SAAO for DoPhot preprocessing.  Added FWHMSKY to
 *                  procfile commands, and GAIN= keyword [rwp/osu]
 * 1.0 1998 May 27: first non-beta release, with documentation [rwp/osu]
 * 1.1 1998 Jun 2:  fixed bug in FWHMSKY header output: 
 *                  TMSAX -> TSMAX [rwp/osu]
 * 1.2 1998 Nov 23: fixed bug in gparams (overwriting array bounds on
 *                  obliterating stars on the working image copy).
 *                  compiled successfully with CFITSIO 2.024 [rwp/osu]
 * 1.3 1999 May 27: fixed more bugs in gparams (sky estimates were 
 *                  noisy and systmatically under-estimated by the
 *                  histogram bin width of 5.33 ADU.  Also tested
 *                  compilation with CFITSIO 2.031 [rwp/osu]
 *  p1 1999 Jun 1:  Weird bug for SHORT, might be due to a deep
 *                  CFITSIO bug (we didn't see it in v1.4), during
 *                  debug I added an argument to FITSError, and 
 *                  the open/close voodoo below.  Weird.  [rwp/osu]
 *
 *  p2 1999 Jul 02: Fixed bugs that caused bogus gain in FWHMSKY (bad
 *                  function prototype), and cleaned up bad history
 *                  cards due to long strings. [rwp/osu]
 *
 * 1.4 1999 Aug 8:  Added DARK feature, in part to support the AP7 CCD
 *                  camera at Perth, but it appears that DARKs are used
 *                  at U Tas as well. [rwp/osu]
 *
 * 1.6.0 2005 Aug 11: added CLIP command to clip data values (after FLAT)
 *                    [rwp/osu]
 *
 ***************************************************************************/

#include "CCDProc.h"

#define VERSION "ccdproc v1.6.0 (2005-08-11)"

struct st st;            /* System table structure */
struct st *systab=&st;   /* set pointers to use the tables globally */

int 
main(int argc, char *argv[])
{
  
  fitsfile *infptr;     /* input image FITS file pointer */ 
  fitsfile *zerofptr;   /* zero image FITS file pointer */
  fitsfile *flatfptr;   /* flat-field image FITS file pointer */
  fitsfile *darkfptr;   /* dark frame FITS file pointer */
  fitsfile *outfptr;    /* output image FITS file pointer */
  
  char hist[10][FLEN_VALUE];     /* history card array */
  int nhist = 0;
  char datetime[FLEN_VALUE];     /* Formatted Date/Time Tag */
  int c;
  int status = 0;
  int nfound, anynull;
  long inaxis = 2;   /* all images are 2-D FITS */
  long znaxis = 2;
  long fnaxis = 2;
  long onaxis = 2;
  long dnaxis = 2;
  long inaxes[2];    /* image dimensions */
  long znaxes[2];
  long fnaxes[2];
  long onaxes[2];
  long dnaxes[2];
  long fpixel = 1;
  long npixels, outpix, i, j;
  float *indata, *zerodata, *flatdata, *outdata, *darkdata;
  float nullval = 0.0; 

  long sc, ec, nc, sr, er, nr;

  char infile[MED_STR_SIZE] = "-";             /* default input is STDIN */
  char outfile[MED_STR_SIZE] = "-";            /* default output is STDOUT */
  char procfile[MED_STR_SIZE] = "ccdproc.pro"; /* default proc file */

  char logstr[MED_STR_SIZE];                /* logfile entry */

  double bzero = 0.0, bscale = 1.0;  /* default scaling factors  */
  double mbias, sigbias;             /* mean and sigma overscan bias */

  float fwhm, sky, tsmin, tsmax;     /* FWHMSKY parameters */

  char errmsg[MED_STR_SIZE];           /* generic error string */

  /* hardcore version: 4 arguments or nothing (incl. command name)*/

  if (argc != 4) {
      Usage(); 
      exit(1); 
    }

  /* parse the command line */

  while (--argc > 0 && (*++argv)[0]=='-')
    {
      c = *++argv[0];
      switch(c) {
      case 'P':
        sscanf(*argv,"P%s",procfile);
        break;
      default:
	sprintf(errmsg,"ccdproc: illegal command-line option %c",c);
	PrintError(errmsg);
        Usage();
        exit(1);
        break;
      }
    }
  if (argc == 1) 
    {
      strcpy(infile,argv[0]);
    } 
  else if (argc == 2)
    {
      strcpy(infile,argv[0]);
      strcpy(outfile,argv[1]);
    }

  /* Open and parse the proc file */

  ParseProc(procfile);

  /* Put in the indentifier in the log */

  if (systab->dolog) {
    BZero(logstr,sizeof(logstr));
    TimeTag(datetime);
    sprintf(logstr,"\n%s: Processing %s using procfile %s\n",datetime,infile,procfile);
    status = write(systab->fd_logfile,logstr,strlen(logstr));
    BZero(logstr,sizeof(logstr));
    sprintf(logstr,"%s: Output filename: %s\n",datetime,outfile);
    status = write(systab->fd_logfile,logstr,strlen(logstr));
    status = 0;
  }

  /********************************************************************/
  /* Open the input FITS file READONLY                                */

  if ( fits_open_file(&infptr, infile, READONLY, &status) )
       FITSError(status,"Cannot open input file");

  /* read the NAXIS1 and NAXIS2 keyword to get image size */

  if ( fits_read_keys_lng(infptr, "NAXIS", 1, 2, inaxes, &nfound, &status) )
       FITSError(status,"Cannot read input file header");

  npixels  = inaxes[0] * inaxes[1];      /* number of pixels in the image */

  /* Allocate memory for the data as a 1-D floating array. */

  indata = (float*) fvalloc(npixels);  

  /* Read in the input file */

  if (fits_read_img(infptr, TFLOAT, fpixel, npixels, &nullval,
                    indata, &anynull, &status))
    FITSError(status,"Cannot read input data");

  /***** Voodoo for SHORT, close then open file again [rwp/osu 1999Jun01 *****/

  if (fits_close_file (infptr, &status))
    FITSError(status,"Cannot close input FITS file");

  if (fits_open_file(&infptr, infile, READONLY, &status) )
       FITSError(status,"Cannot open input file");

  /***** this takes care of stale info from doing read_keys etc. for shorts *****/

  /* if we are doing overscan bias correction, do it now */

  if (systab->dobias) {

    for (i=0; i<(systab->namps); i++) {

      /* compute the mean bias in the current bias section */

      mbias = immean(indata,inaxes[0],inaxes[1],systab->bsc[i],
		     systab->bec[i],systab->bsr[i],systab->ber[i],&sigbias);

      /* subtract the mean bias from the current data section */

      sc = systab->dsc[i] - 1;
      sr = systab->dsr[i] - 1;
      nc = systab->dec[i] - systab->dsc[i] + 1;
      nr = systab->der[i] - systab->dsr[i] + 1;
	
      AddBox(indata,inaxes[0],sc,nc,sr,nr,-mbias);

      /* make a history card for the final output header */

      TimeTag(datetime);
      sprintf(hist[nhist],"(ccdproc:%s) Overscan %d of %d: %.2f+/-%.2f",
	      datetime,i+1,systab->namps,mbias,sigbias);

      if (systab->dolog) {
	BZero(logstr,sizeof(logstr));
	sprintf(logstr,"%s\n",hist[nhist]);
	status = write(systab->fd_logfile,logstr,strlen(logstr));
	status = 0;
      }

      nhist++;
    }

  } // end (systab->dobias)

  /* 
   * If we are doing trimming, set the limits appropriately, allocate VM
   * for the subarray, and do the trimming.  
   */

  if (systab->dotrim) {

    nc = systab->tec - systab->tsc + 1;
    sc = systab->tsc;
    ec = systab->tec;
    nr = systab->ter - systab->tsr + 1;
    sr = systab->tsr;
    er = systab->ter;
    npixels = nc*nr;
    outdata = (float*) fvalloc(npixels);  

    if (imtrim(indata,inaxes[0],inaxes[1],outdata,sc,ec,sr,er)) {
      PrintError("Error: ccdproc could not trim image");
      exit(SYSERR);
    }

    /* make a history card for the final output header */
    
    TimeTag(datetime);
    sprintf(hist[nhist],"(ccdproc:%s) Trimmed to [%d:%d,%d:%d]",
	    datetime,sc,ec,sr,er);

    if (systab->dolog) {
      BZero(logstr,sizeof(logstr));
      sprintf(logstr,"%s\n",hist[nhist]);
      status = write(systab->fd_logfile,logstr,strlen(logstr));
      status = 0;
    }

    nhist++;

    /* 
     * discard the input data array as we no longer need it, but leave the
     * FITS file open so we can copy the header into the output file below.
     */

    free(indata);

  } 
  else {  // use the entire image, untrimmed

    nc = inaxes[0];
    sc = 1;
    ec = nc;
    nr = inaxes[1];
    sr = 1;
    er = nr;
    npixels = nc*nr;

  } // end systab->dotrim

  /********************************************************************
   * If doing zero subtraction, Open the zero file, sniff the header for
   * the size and die if not the same size as the data array.
   * Subtract the zero frame from the data array, then close the Zero
   * file and free its virtual memory.
   */

  if (systab->dozero) {

    if ( fits_open_file(&zerofptr, systab->zerofile, READONLY, &status) )
      FITSError(status,"Cannot open zero frame");

    /* read the NAXIS1 and NAXIS2 keyword to get image size */

    if ( fits_read_keys_lng(zerofptr,"NAXIS",1,2,znaxes,&nfound,&status) )
      FITSError(status,"Cannot read zero frame header");

    if ((znaxes[0] != nc) || (znaxes[1] != nr)) {
      PrintError("Error: Zero and Data Images are not the same size");
      exit(SYSERR); 
    }

    /* It's the same size, allocate a 1-D floating array for its pixels */

    zerodata = (float*) fvalloc(npixels);  

    /* Read in the zero file */

    if (fits_read_img(zerofptr, TFLOAT, fpixel, npixels, &nullval,
		      zerodata, &anynull, &status))
      FITSError(status,"Cannot read zero frame data");

    /* Subtract the zero frame from the image (trimmed or original) */

    if (systab->dotrim)
      ImSub(outdata, zerodata, npixels);
    else
      ImSub(indata, zerodata, npixels);

    /* Close the zero file & free its VM */

    if (fits_close_file (zerofptr, &status))
      FITSError(status,"Cannot close zero frame");

    free(zerodata);

    /* make a history card for the final output header */
    
    TimeTag(datetime);
    sprintf(hist[nhist],"(ccdproc:%s) Zero: %s",
	    datetime,systab->zerofile);

    if (systab->dolog) {
      BZero(logstr,sizeof(logstr));
      sprintf(logstr,"%s\n",hist[nhist]);
      status = write(systab->fd_logfile,logstr,strlen(logstr));
      status = 0;
    }

    nhist++;

  } // end systab->dozero

  /********************************************************************
   * If doing dark subtraction, Open the dark file, sniff the header for
   * the size and die if not the same size as the data array.
   * Subtract the dark frame from the data array, then close the dark
   * file and free its virtual memory.  
   *
   * We are implicitly assuming that the dark frame has been properly
   * pre-processed, at least OT or OTZ as required.  In this reduction
   * model dark subtraction comes before flat-field correction, but
   * after both overscan bias/trim and zero (2-D bias) subtraction.
   */

  if (systab->dodark) {

    if ( fits_open_file(&darkfptr, systab->darkfile, READONLY, &status) )
      FITSError(status,"Cannot open dark frame");

    /* read the NAXIS1 and NAXIS2 keyword to get image size */

    if ( fits_read_keys_lng(darkfptr,"NAXIS",1,2,dnaxes,&nfound,&status) )
      FITSError(status,"Cannot read dark frame header");

    if ((dnaxes[0] != nc) || (dnaxes[1] != nr)) {
      PrintError("Error: Dark and Data Images are not the same size");
      exit(SYSERR); 
    }

    /* It's the same size, allocate a 1-D floating array for its pixels */

    darkdata = (float*) fvalloc(npixels);  

    /* Read in the dark file */

    if (fits_read_img(darkfptr, TFLOAT, fpixel, npixels, &nullval,
		      darkdata, &anynull, &status))
      FITSError(status,"Cannot read dark frame data");

    /* Scale the dark frame if systab->darkscale is not 1.0 */

    if (systab->darkscale != 1.0) {
      MulCon(darkdata, systab->darkscale, npixels);
    }

    /* Subtract the dark frame from the image */

    if (systab->dotrim)
      ImSub(outdata, darkdata, npixels);
    else
      ImSub(indata, darkdata, npixels);

    /* Close the dark file & free its VM */

    if (fits_close_file (darkfptr, &status))
      FITSError(status,"Cannot close dark frame");

    free(darkdata);

    /* make a history card for the final output header */
    
    TimeTag(datetime);
    sprintf(hist[nhist],"(ccdproc:%s) Dark: %s scale=%.2f",
	    datetime,systab->darkfile,systab->darkscale);

    if (systab->dolog) {
      BZero(logstr,sizeof(logstr));
      sprintf(logstr,"%s\n",hist[nhist]);
      status = write(systab->fd_logfile,logstr,strlen(logstr));
      status = 0;
    }

    nhist++;

  } // end systab->dodark

  /********************************************************************
   * If doing flat-field division, open the flat file, sniff the header for
   * the size and die if not the same size as the data array.
   * Divide the data array by the flat field.  We assume the flat has
   * been properly normalized.  Once finished, close the flat file and
   * free its virtual memory.
   */

  if (systab->doflat) {
        
    if ( fits_open_file(&flatfptr, systab->flatfile, READONLY, &status) )
      FITSError(status,"Cannot open flat field frame");

    /* read the NAXIS1 and NAXIS2 keyword to get image size */

    if ( fits_read_keys_lng(flatfptr,"NAXIS",1,2,fnaxes,&nfound,&status) )
      FITSError(status,"Cannot read flat field frame header");

    if ((fnaxes[0] != nc) || (fnaxes[1] != nr)) {
      PrintError("Error: Flat Field and Data Images are not the same size");
      exit(SYSERR); 
    }

    /* It's the same size, allocate a 1-D floating array for its pixels */

    flatdata = (float*) fvalloc(npixels);  

    /* Read in the flat file */

    if (fits_read_img(flatfptr, TFLOAT, fpixel, npixels, &nullval,
		      flatdata, &anynull, &status))
      FITSError(status,"Cannot read flat field frame data");

    /* Divide by the flat field */

    if (systab->dotrim)
      ImDiv(outdata, flatdata, npixels);
    else
      ImDiv(indata, flatdata, npixels);

    /* Close the flat file & free its VM */

    if (fits_close_file (flatfptr, &status))
      FITSError(status,"Cannot close flat field frame");

    free(flatdata);

    /* make a history card for the final output header */
    
    TimeTag(datetime);
    sprintf(hist[nhist],"(ccdproc:%s) Flat: %s",
	    datetime,systab->flatfile);

    if (systab->dolog) {
      BZero(logstr,sizeof(logstr));
      sprintf(logstr,"%s\n",hist[nhist]);
      status = write(systab->fd_logfile,logstr,strlen(logstr));
      status = 0;
    }

    nhist++;

  } // end systab->doflat



  // If doing clipping, do it here, before any further computational
  // steps (e.g., FWHM/Sky or scaling for output data conversion)

  if (systab->doclip) {
    if (systab->dotrim)
      ImClip(outdata, systab->minclip, systab->maxclip, npixels);
    else
      ImClip(indata, systab->minclip, systab->maxclip, npixels);

    /* make a history card for the final output header */
    
    TimeTag(datetime);
    sprintf(hist[nhist],"(ccdproc:%s) CLIP between %.2f..%.2f",
	    datetime, systab->minclip, systab->maxclip);

    if (systab->dolog) {
      BZero(logstr,sizeof(logstr));
      sprintf(logstr,"%s\n",hist[nhist]);
      status = write(systab->fd_logfile,logstr,strlen(logstr));
      status = 0;
    }
    nhist++;
  }

  /*****************************************************************
   * Compute FWHM and Sky parameters for DoPhot if required.  These
   * data will be placed in the FITS header.
   */

  if (systab->dofwhm) {
    if (systab->dotrim)
      gparam(nc,nr,outdata,systab->detgain,&fwhm,&sky,&tsmin,&tsmax);
    else
      gparam(nc,nr,indata,systab->detgain,&fwhm,&sky,&tsmin,&tsmax);

    /* If fwhm = -1, it means no stars were found */

    if (fwhm == -1.0) 
      PrintError("Warning: FWHMSKY found no stars");

    /* make a history card for the final output header */
    
    TimeTag(datetime);
    sprintf(hist[nhist],"(ccdproc:%s) FHWMSKY with GAIN=%.2f",
	    datetime,systab->detgain);

    if (systab->dolog) {
      BZero(logstr,sizeof(logstr));
      sprintf(logstr,"%s\n",hist[nhist]);
      status = write(systab->fd_logfile,logstr,strlen(logstr));
      status = 0;
    }

    nhist++;

  } // end systab->dofwhm

  /*****************************************************************
   * Prepare the output FITS file 
   *
   * compute the scaling parameters (BZERO & BSCALE) if creating
   * an integer FITS file.  These are unnecesary for floats
   */

  if (systab->bitpix != FLOAT_IMG) {
    if (systab->dotrim) {
      if (fitscl(outdata,npixels,systab->bitpix,&bzero,&bscale)) {
	PrintError("Error: ccdproc could not compute BZERO/BSCALE");
	exit(SYSERR);
      }
    } else {
      if (fitscl(indata,npixels,systab->bitpix,&bzero,&bscale)) {
	PrintError("Error: ccdproc could not compute BZERO/BSCALE");
	exit(SYSERR);
      }
    }
  }

  /* Open a FITS file for the output */

  if (fits_create_file (&outfptr, outfile, &status))
    FITSError(status,"Cannot create output FITS file");

  /* copy the header from the input image */

  if (fits_copy_hdu (infptr, outfptr, 32, &status))
    FITSError(status,"Cannot copy header to output FITS file");

  /* Remove any BZERO or BSCALE cards (ignore errors) */

  fits_delete_key (outfptr,"BZERO",&status);
  status = 0;
  fits_delete_key (outfptr,"BSCALE",&status);
  status = 0;

  /* Modify the data type of the current array */

  onaxes[0] = nc;
  onaxes[1] = nr;

  if (fits_resize_img(outfptr,systab->bitpix,onaxis,onaxes,&status))
    FITSError(status,"Cannot convert data for output FITS file");

  /* If doing integers, set the BZERO and BSCALE parameters */

  if (systab->bitpix != FLOAT_IMG) {
    if (fits_set_bscale(outfptr,bscale, bzero, &status))
      FITSError(status,"Cannot set output FITS header scaling params"); 

    /* write the BZERO and BSCALE cards into the header (why CFITSIO doesn't
       do this when it scales we don't know) */

    if (fits_write_key_dbl(outfptr,"BSCALE",bscale,14,
			   "DATA=FITS*BSCALE+BZERO",&status))
      FITSError(status,"Cannot write BSCALE card"); 

    if (fits_write_key_dbl(outfptr,"BZERO",bzero,14,
			   "see BSCALE",&status))
      FITSError(status,"Cannot write BZERO card"); 

  }

  /* write the data to the output image, datatype conversion is implicit */

  status = 0;
  if (systab->dotrim) {
    fits_write_img (outfptr, TFLOAT, 1, npixels, outdata, &status);
    if (status !=0) 
      FITSError(status,"Cannot write trimmed output FITS image"); 
  } else {
    fits_write_img (outfptr, TFLOAT, 1, npixels, indata, &status);
    if (status !=0) 
      FITSError(status,"Cannot write output FITS image"); 
  }

  /* If we computed FWHM & SKY parameters, put them in the header */

  if (systab->dofwhm) {
    if (fits_write_key_flt(outfptr,"SKY",sky,10,
			   "Modal SKY level (ADU)",&status))
      FITSError(status,"Cannot write SKY card"); 

    if (fits_write_key_flt(outfptr,"FWHM",fwhm,10,
			   "Mean stellar FWHM (pix)",&status))
      FITSError(status,"Cannot write FWHM card"); 

    if (fits_write_key_flt(outfptr,"TSMIN",tsmin,10,
			   "DoPhot min threshold (THRESHMIN)",&status))
      FITSError(status,"Cannot write TSMIN card"); 

    if (fits_write_key_flt(outfptr,"TSMAX",tsmax,10,
			   "DoPhot max threshold (THRESHMAX)",&status))
      FITSError(status,"Cannot write TSMAX card"); 
  }

  /* add history cards to the FITS header recording what was done */

  if (nhist > 0) {
    for (i=0;i<nhist;i++) {
      if (fits_write_history (outfptr, hist[i], &status))
	FITSError(status,"Cannot write HISTORY card(s)");
    }
  }

  /* update the runtime log and then close it */

  if (systab->dolog) {
    TimeTag(datetime);
    BZero(logstr,sizeof(logstr));
    sprintf(logstr,"(ccdproc:%s) Done\n",datetime);
    status = write(systab->fd_logfile,logstr,strlen(logstr));
    status = close(systab->fd_logfile);
  }

  /* Done! close the remaining open FITS files */

  status = 0;

  if (fits_close_file (infptr, &status))
    FITSError(status,"Cannot close input FITS file");

  if (fits_close_file (outfptr, &status))
    FITSError(status,"Cannot close output FITS file");
  
  return 0;

}

/*************************************************************************** 
 *
 *  PrintError - Print generic error messages w/o exiting.
 *
 */ 

void PrintError(char *errstr)
{
  char logstr[MED_STR_SIZE];
  int ierr = 0;

  fprintf(stderr,"%s\n",errstr);
  if (systab->dolog) {
    BZero(logstr,sizeof(logstr));
    sprintf(logstr,"   %s\n",errstr);
    ierr = write(systab->fd_logfile,logstr,strlen(logstr));
  }    

}

/***************************************************************************/

void Usage()
{
  fprintf(stderr, "\nUsage: ccdproc -Pprocfile inFITS outFITS\n");
  fprintf(stderr, "  where procfile = the process control file\n");
  fprintf(stderr, "        inFITS  = input FITS filename\n");
  fprintf(stderr, "        outFITS = output FITS filefile\n");
  fprintf(stderr, "see ccdproc.hlp for details\n");
}

/*
 *  Allocate memory for a 1-dimensional array of floats
 */

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
 * ImSub - subtract two images 
 *
 * arguments:
 *   im1 = image 1 (float)
 *   im2 = image 2 (float)
 *     n = total number of pixels in both images (long)
 *
 * returns:
 *   replaces im1 by im1-im2
 * 
 * notes:
 *   images must be perfectly aligned (same number of pixels and formats)
 *
 ***************************************************************************/

void ImSub (float *im1, float *im2, long n)
{
  int i;
  register float *a, *b; 
  long N;
  
  a = im1; 
  b = im2; 
  N = n; 

  for (i=0; i<N; i++)
    *a++ -= *b++; 

}

/***************************************************************************
 *
 * ImDiv - divide two images with zero supression
 *
 * arguments:
 *   im1 = image 1 (float)
 *   im2 = image 2 (float)
 *     n = total number of pixels in both images (long)
 *
 * returns:
 *   replaces im1 by im1/im2, division by 0 results in 0.
 * 
 * notes:
 *   images must be perfectly aligned (same number of pixels and formats)
 *
 ***************************************************************************/

void ImDiv (float *im1, float *im2, long n)
{
  int i;
  register float *a, *b; 
  long N;

  a = im1; 
  b = im2; 
  N = n; 

  for (i=0; i<N; i++) {
    if (*b == 0.0 ) 
      *a++ = 0.0; 
    else
      *a++ /= *b;
    *b++;
  }

}


/***************************************************************************
 *
 * ImClip - Clip image data values outside a given range
 *
 * arguments:
 *   img = image (float)
 *   minval = minimum data value (float)
 *   maxval = maximum data value (float)
 *     n = number of pixels (long)
 *
 * If a pixel has data<min, replace it with min, if data>max, replaces it
 * with max, otherwise it remains unchanged
 * 
 ***************************************************************************/

void 
ImClip(float *img, float minval, float maxval, long n)
{
  int i;
  float *a;
  
  a = img; 

  for (i=0; i<n; i++) {
    if (*a < minval)
      *a = minval;
    else if (*a > maxval)
      *a = maxval;
    *a++;
  }
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
 * AddBox - add constant to image subarray
 *
 * Arguments:
 *    a (register/float) = array with data (full array)
 *    ncola (int) = number of columns in full array
 *    sc, ncol (int) = starting col & ncols of subarray
 *    sr, nrow (int) = starting row & nrows of subarray
 *    f (float) = constant to add
 *
 * Adapted from addcon by Tod Lauer, see also AddCon()
 *
 ***************************************************************************/

void AddBox (float *a, int ncola, int sc, int ncol, int sr, int nrow,
	     float f)
{
  register float *aim, c;
  int irow, icol;

  aim     = a + sr*ncola + sc;    /* Initialize pointer   */
  c       = f;
  for (irow = 1; irow <= nrow; irow++) {
    for (icol = 1; icol <= ncol; icol++)
      *aim++  += c;
    aim += ncola - ncol;
  }
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
 * Modified 1999 Jul 2 to generate ISO8604-compliant Y2K format tags.
 *
 ***************************************************************************/

void 
TimeTag(char *tagstr)
{
  time_t now;
  struct tm *tm;

  now = time(NULL);
  tm = localtime(&now);
  /*  if (strftime(tagstr,FLEN_VALUE,"%c",tm)==0) */
  if (strftime(tagstr,FLEN_VALUE,"%Y-%m-%dT%H:%M:%S",tm)==0)
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
 * 16 and 32-bit precision min/max are hardwired to 1-2**(BITPIX-1) and
 * 2**(BITPIX-1)-1, respectively.
 *
 * R. Pogge & P. Martini
 * OSU Astronomy
 * 1998 April 28
 *
 ***************************************************************************/

int fitscl(float *array, long npix, int bitpix, double *bzero, double *bscale)
{

  long ii;
  double dmin, dmax;
  double pixel;
  double fitsmin, fitsmax;

  /* get data limits for bitpix.  If bitpix is not 16 or 32, error */

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

  /* compute the data min and max values */

  dmin  = 1.0E30;
  dmax  = -1.0E30;

  for (ii = 0; ii < npix; ii++)  {
    pixel = (double)(array[ii]);
    if ( pixel < dmin )
      dmin = pixel;
    
    if ( pixel > dmax )
      dmax = pixel;
  }

  /* compute BZERO and BSCALE appropriate for BITPIX */

  *bzero = ((dmin*fitsmax)-(dmax*fitsmin))/(fitsmax-fitsmin);
  if (dmax == dmin) {
    *bscale = 1.0;
  } else {
    *bscale = (dmax-dmin) / (fitsmax-fitsmin);
  }
  return(0);
  
}


/***************************************************************************
 *
 * ParseProc - Parse the process control (proc) file and set up the
 *             keyword tables
 * 
 * Requirements: none
 *
 * Returns: none
 *
 * R. Pogge
 * OSU Astronomy Dept.
 * 1998 May 19
 *
 * Based on the ini file parser from caliban by B. Hartung
 * See ccdproc.hlp for the list of keywords and their syntax
 *
 ***************************************************************************/

void ParseProc(char *procfile)
{
  int lcv;                             /* Loop control variable */
  int ii;                              /* generic integer */
  char keyword[MED_STR_SIZE];          /* Header Keywords  */
  char argbuf[MED_STR_SIZE];           /* Generic argument buffer */
  char workbuf[MED_STR_SIZE];          /* Generic working buffer */
  char inbuf[BUF_SIZE];                /* Generic input buffer */
  char logfile[MED_STR_SIZE];          /* runtime logfile name string */
  FILE *rfp;                           /* Initialization file pointer */

  struct tm *loctm;
  time_t now;

  /* Some initializations */

  systab->namps = 0;
  systab->detgain = 1.0;
  systab->dobias = ccd_FALSE;
  systab->dozero = ccd_FALSE;
  systab->dodark = ccd_FALSE;
  systab->doflat = ccd_FALSE;
  systab->dotrim = ccd_FALSE;
  systab->dolog  = ccd_FALSE;
  systab->bitpix = FLOAT_IMG;
  systab->dofwhm = ccd_FALSE;
  systab->darkscale = 1.0;

  /* Open the initialization file, exit abnormally if absent */

  if(!(rfp=fopen(procfile, "r"))) 
    {
      fprintf(stderr,"Error: ccdproc unable to open procfile %s\n",procfile);
      fprintf(stderr,"Reason: %s\n", ERRORSTR);
      exit(SYSERR);
    }

  /* 
   * Loop through file reading 80 byte blocks at a time.  Skip 
   * Skip comments (#'s) and blank lines as you go.
   *
   * Note that not all possible keywords are relevant for all
   * programs.  We do this to provide a common syntax to avoid
   * unnecessary replication of this block of code.
   */

  else
    {

      while(fgets(inbuf, 80, rfp))
	{
	  if((inbuf[0]!='#') && (inbuf[0]!='\n') && inbuf[0]!=NUL) 
	    {
	      inbuf[80] = NUL;

	      /* 
	       * Convert to keyword uppercase to ease the comparison, values
	       * of keywords are taken without case conversion
	       */

	      GetArg(inbuf, 1, argbuf);
	      strcpy(keyword, argbuf);
	      UpperCase(keyword);
	      
	      /* DTYPE: output data type (SHORT|LONG|FLOAT[default]) */

	      if (strcmp(keyword, "DTYPE")==0)
		{
		  GetArg(inbuf, 2, argbuf);
		  strcpy(workbuf,argbuf);
		  UpperCase(workbuf);
		  if (strcmp(workbuf,"SHORT")==0) 
		    {
		      systab->bitpix = SHORT_IMG;
		    }
		  else if (strcmp(workbuf,"LONG")==0)
		    {
		      systab->bitpix = LONG_IMG;
		    }
		  else
		    {
		      systab->bitpix = FLOAT_IMG;
		    }
		}

	      else if (strcmp(keyword, "LOGFILE")==0)
		{
		  GetArg(inbuf, 2, argbuf);
		  strcpy(workbuf,argbuf);
		  UpperCase(workbuf);

		  /* DATE = ccdproc.yymmdd.log, using the local time */

		  if (strcmp(workbuf,"DATE")==0)
		    {
		      now = time(NULL);
		      loctm = localtime(&now);
		      strftime(logfile,MED_STR_SIZE,"ccdproc.%y%m%d.log",loctm);
		      systab->dolog = ccd_TRUE;
		      systab->fd_logfile=open(logfile,(O_WRONLY|O_CREAT|O_APPEND));
		      if (systab->fd_logfile == -1) {
			fprintf(stderr,"Error: ccdproc cannot open logfile %s\n",logfile);
			fprintf(stderr,"Reason: %s\n", ERRORSTR);
			exit(SYSERR);
		      }
		      chmod(logfile,0666);
		    }

		  /* UTDATE = ccdproc.yymmdd.log, using the UT time */

		  else if (strcmp(workbuf,"UTDATE")==0)
		    {
		      now = time(NULL);
		      loctm = gmtime(&now);
		      strftime(logfile,MED_STR_SIZE,"ccdproc.%y%m%d.log",loctm);
		      systab->dolog = ccd_TRUE;
		      systab->fd_logfile=open(logfile,(O_WRONLY|O_CREAT|O_APPEND));
		      if (systab->fd_logfile == -1) {
			fprintf(stderr,"Error: ccdproc cannot open logfile %s\n",logfile);
			fprintf(stderr,"Reason: %s\n", ERRORSTR);
			exit(SYSERR);
		      }
		      chmod(logfile,0666);
		    }

		  /* Use the argument as the logfile name */

		  else
		    {
		      systab->dolog = ccd_TRUE;
		      systab->fd_logfile=open(argbuf,(O_WRONLY|O_CREAT|O_APPEND));
		      if (systab->fd_logfile == -1) {
			fprintf(stderr,"Error: ccproc cannot open logfile %s\n",logfile);
			fprintf(stderr,"Reason: %s\n", ERRORSTR);
			exit(SYSERR);
		      }
		      chmod(argbuf,0666);
		    }
		}

	      /* BIASSEC: coordinates of bias section (default: none) */
	      /* BIASSEC1: bias section of amp 1 */

	      else if (strcmp(keyword,"BIASSEC")==0 || 
		       strcmp(keyword,"BIASSEC1")==0 )
		{
		  systab->dobias = ccd_TRUE;
		  systab->namps++;
		  GetArg(inbuf,2,argbuf);
		  sscanf(argbuf,"[%d:%d,%d:%d]",&systab->bsc[0],
			 &systab->bec[0],&systab->bsr[0],&systab->ber[0]);
		}

	      /* BIASSEC2: bias section of amp 2 */

	      else if (strcmp(keyword,"BIASSEC2")==0 )
		{
		  systab->dobias = ccd_TRUE;
		  systab->namps++;
		  GetArg(inbuf,2,argbuf);
		  sscanf(argbuf,"[%d:%d,%d:%d]",&systab->bsc[1],
			 &systab->bec[1],&systab->bsr[1],&systab->ber[1]);
		}

	      /* BIASSEC3: bias section of amp 3 */

	      else if (strcmp(keyword,"BIASSEC3")==0 )
		{
		  systab->dobias = ccd_TRUE;
		  systab->namps++;
		  GetArg(inbuf,2,argbuf);
		  sscanf(argbuf,"[%d:%d,%d:%d]",&systab->bsc[2],
			 &systab->bec[2],&systab->bsr[2],&systab->ber[2]);
		}

	      /* BIASSEC4: bias section of amp 4 */

	      else if (strcmp(keyword,"BIASSEC4")==0 )
		{
		  systab->dobias = ccd_TRUE;
		  systab->namps++;
		  GetArg(inbuf,2,argbuf);
		  sscanf(argbuf,"[%d:%d,%d:%d]",&systab->bsc[3],
			 &systab->bec[3],&systab->bsr[3],&systab->ber[3]);
		}

	      /* DATASEC: coordinates of Data section (default: none) */
	      /* DATASEC1: Data section of amp 1 */

	      else if (strcmp(keyword,"DATASEC")==0 || 
		       strcmp(keyword,"DATASEC1")==0 )
		{
		  GetArg(inbuf,2,argbuf);
		  sscanf(argbuf,"[%d:%d,%d:%d]",&systab->dsc[0],
			 &systab->dec[0],&systab->dsr[0],&systab->der[0]);
		}

	      /* DATASEC2: Data section of amp 2 */

	      else if (strcmp(keyword,"DATASEC2")==0)
		{
		  GetArg(inbuf,2,argbuf);
		  sscanf(argbuf,"[%d:%d,%d:%d]",&systab->dsc[1],
			 &systab->dec[1],&systab->dsr[1],&systab->der[1]);
		}

	      /* DATASEC3: Data section of amp 3 */

	      else if (strcmp(keyword,"DATASEC3")==0)
		{
		  GetArg(inbuf,2,argbuf);
		  sscanf(argbuf,"[%d:%d,%d:%d]",&systab->dsc[2],
			 &systab->dec[2],&systab->dsr[2],&systab->der[2]);
		}

	      /* DATASEC4: Data section of amp 4 */

	      else if (strcmp(keyword,"DATASEC4")==0)
		{
		  GetArg(inbuf,2,argbuf);
		  sscanf(argbuf,"[%d:%d,%d:%d]",&systab->dsc[3],
			 &systab->dec[3],&systab->dsr[3],&systab->der[3]);
		}

	      /* ZERO: Zero (2-D bias) FITS file */

	      else if (strcmp(keyword, "ZERO")==0) 
		{
		  systab->dozero = ccd_TRUE;
		  GetArg(inbuf, 2, argbuf);
		  strcpy(systab->zerofile,argbuf);
		}

	      /* DARK: Dark frame FITS file and the scaling factor.
	       *       The default scaling factor is 1.0 if absent 
               */
	      else if (strcmp(keyword, "DARK")==0) 
		{
		  systab->dodark = ccd_TRUE;
		  GetArg(inbuf, 2, argbuf);
		  strcpy(systab->darkfile,argbuf);
		  GetArg(inbuf, 3, argbuf);
		  if (argbuf[0] == NUL) {
		    systab->darkscale = 1.0 ;
		  } else {
		    sscanf(argbuf,"%f",&systab->darkscale);
		  }
		}

	      /* FLAT: Flat-field FITS file */

	      else if (strcmp(keyword, "FLAT")==0) 
		{
		  systab->doflat = ccd_TRUE;
		  GetArg(inbuf, 2, argbuf);
		  strcpy(systab->flatfile,argbuf);
		}

	      // CLIP: Clip image values outside a given range (default: none)
	      //       Syntax: CLIP min max
	      //       Example: CLIP 0 65000

	      else if (strcmp(keyword,"CLIP")==0) {
		systab->doclip = ccd_TRUE; 
		GetArg(inbuf,2,argbuf);
		systab->minclip = (float)(atof(argbuf));
		GetArg(inbuf,3,argbuf);
		systab->maxclip = (float)(atof(argbuf));
	      }

	      /* TRIM: coordinates of trim section (default: none) */

	      else if (strcmp(keyword,"TRIM")==0)
		{
		  systab->dotrim = ccd_TRUE;
		  GetArg(inbuf,2,argbuf);
		  sscanf(argbuf,"[%d:%d,%d:%d]",&systab->tsc,&systab->tec,
			 &systab->tsr,&systab->ter);
		}

	      /* 
	       * FWHMSKY: compute FWHM of stars & SKY level in the reduced image
	       *       Optional second argument GAIN= gives the detector gain
	       *       in electrons/ADU (default: 2 e-/ADU) 
               */

	      else if (strcmp(keyword,"FWHMSKY")==0)
		{
		  systab->dofwhm = ccd_TRUE;
		  GetArg(inbuf,2,argbuf);
		  if (argbuf[0] != '0') {
		    strcpy(workbuf,argbuf);
		    UpperCase(workbuf);
		    sscanf(workbuf,"GAIN=%f",&systab->detgain);
		  }
		}

	      /* unrecognized directive in the proc file, squawk & exit */

	      else 
		{
		  fprintf(stderr,"Error: Unrecognized keyword %s in procfile %s\n",
			 keyword,procfile);
		  exit(SYSERR);
		}
	    }
	  BZero(inbuf, sizeof(inbuf)); /* Reset input buffer */
	}
    }
  
  if (rfp!=0)
    fclose(rfp);              /* Close initialization file */

}

/***************************************************************************
 *
 * imtrim - trim an array down to a given size
 *
 * Arguments:
 *   inarray (float): original data array of size nc*nr (1-D array)
 *   nc, nr (long): dimensions of the input array
 *   outarray (float): output subarray (1-D array of correct subarray size)
 *   sc,ec (long): starting & ending column pixels for the subarray (axis1)
 *   sr,er (long): starting & ending row pixels for the subarray (axis2)
 *
 * R. Pogge
 * OSU Astronomy
 * 1998 May 18
 *
 ***************************************************************************/

int imtrim(float *inarray, long nc, long nr, float *outarray, 
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

/***************************************************************************
 *
 * immean - compute mean and sigma in an image or subarray
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

  if (systab->dolog) {
    BZero(logstr,sizeof(logstr));
    sprintf(logstr,"   ccdproc terminated on fatal cfitsio error\n");
    ierr = write(systab->fd_logfile,logstr,strlen(logstr));
  }    

  fits_get_errstatus(status, status_str);   /* get the error description */
  fprintf(stderr, "status=%d: %s\n", status, status_str);

  if (systab->dolog) {
    BZero(logstr,sizeof(logstr));
    sprintf(logstr,"   Error %d: %s\n",status,status_str);
    ierr = write(systab->fd_logfile,logstr,strlen(logstr));
  }    

  /* get first message; null if stack is empty */
  if ( fits_read_errmsg(fitserrmsg) )
    {
      fprintf(stderr, " %s\n", fitserrmsg);
      if (systab->dolog) {
	BZero(logstr,sizeof(logstr));
	sprintf(logstr,"   %s\n",fitserrmsg);
	ierr = write(systab->fd_logfile,logstr,strlen(logstr));
      }    
      
      while ( fits_read_errmsg(fitserrmsg) )  { /* get remaining messages */
        fprintf(stderr, " %s\n", fitserrmsg);
	if (systab->dolog) {
	  BZero(logstr,sizeof(logstr));
	  sprintf(logstr,"   %s\n",fitserrmsg);
	  ierr = write(systab->fd_logfile,logstr,strlen(logstr));
	}    
      }
    }

  exit(status);       /* terminate the program, returning error status */

}

