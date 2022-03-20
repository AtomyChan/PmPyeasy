// 
// fwhmsky -- Compute DoPhot fwhm, sky, and threshold values.
//
// Include file for fwhmsky.c
//
// R. Pogge
// OSU Astronomy Dept.
// 1999 May 
//**************************************************************************/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "fitsio.h"    // CFITSIO package header file

//
// Select the detector system, as this defines the parameters used
// by FWHMsky - we need to do this better someday...
//

#define Fairchild477
#undef  KELT4K

// Detector Gain in electrons/ADU 
// gain is 8 e-/ADU for ANDICAM before 2000 June 1, based on measurements 
// gain is 2.3 e-/ADU for the new Fairchild 477 CCD [2003 Feb]
//
// The maximum data value to be used for this type of CCD image
//
// New ANDICAM Fairchild 477 CCD: 60000


#if defined(Fairchild477)
#define DETGAIN 2.3
#define MAXDATA 60000
#elif defined (KELT4K)
#define DETGAIN 3.6
#define MAXDATA 12000
#endif

// Define __DEBUG if you want verbose output

#undef __DEBUG

// Function Prototypes

void FWHMError(int status);
void FWHMUsage(void);
float *fvalloc(int N);

void gparam(int ,int ,float *,float ,float *,float *,float *,float *);
void gparam2(int ,int ,float *,float ,float *,float *,float *,float *,float *);

void getsky(int ,int ,float *,int *,float *);

void getfwhm(int ,int ,int ,float ,int *,float *,float *, float *, float );
void getfwhm2(int ,int ,int ,float ,int *,float *,float *, float *, float *, float );

void fit2(int ,int ,int ,int ,float *,float *,int *);

// maximum number of stars to use for FWHM stats 

#define MAXSTARS 10

// maximum number of saturated pixels 

#define MAXSAT 10000

static int magic = INT_MAX;
static float bot = 1.0; // was -1000.0;
static float top = MAXDATA;
static float pstar = 1.0;
static int cmax = MAXDATA;
static int starmax = MAXDATA;
static int ibor = 30;

// changed nhist to be MAXDATA, defined in fwhmsky.h

static int nhist=MAXDATA;

