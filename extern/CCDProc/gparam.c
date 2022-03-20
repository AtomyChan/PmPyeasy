//
// 
// gparam.c - routines used by fwhmsky.c to estimate starting 
//            values for DoPhot   
//
// Author:
//   John Menzies (SAAO)
//
// Modified by R. Pogge (OSU) for use with the local PLANET pipeline and
// the MicroFUN pipelines (2002 June)
//
//---------------------------------------------------------------------------

#include <stdio.h>
#include <malloc.h>
#include <limits.h>

#include "fwhmsky.h"  // package header file

//---------------------------------------------------------------------------

void 
gparam(int naxis1,int naxis2,float *indata,float gain,float *fwhm,
       float *sky,float *tsmin,float *tsmax)
{
  float sgx[MAXSTARS],sgy[MAXSTARS],sacc[7],slim[7] ;
  int *noise ;
  int i,j,k,jx,jy,kk,i1,i2,j1,j2,mx ;

  noise = (int *)(malloc(naxis1*naxis2*sizeof(int)));

  getsky(naxis1,naxis2,indata,noise,sky) ;
#ifdef __DEBUG
  printf("Sky: %f\n",*sky) ;
#endif
  getfwhm(naxis1,naxis2,ibor,gain,noise,fwhm,tsmin,tsmax,*sky) ;

  free(noise) ;
}

//---------------------------------------------------------------------------
          
void 
getsky(int naxis1,int naxis2,float *indata,int *noise,float *sky)
{
  float sgx[MAXSTARS],sgy[MAXSTARS],sacc[7],slim[7] ;
  float aa,bb,sz,szz ;
  int isat[MAXSAT],jsat[MAXSAT],nsat,npix,z ;
  int *hist ;
  int i,j,k,jx,jy,kk,i1,i2,j1,j2,mx ;

  // allocate space for histogram of counts 

  hist = (int *)(calloc(nhist,sizeof(int))) ;
  aa = (float)(nhist) / (top-bot) ;
  bb = -aa * bot ;

  // fill histogram, between limits bot,top; locate saturated pixels 

  nsat = npix = 0 ;
  isat[nsat] = -naxis1 ;
  jsat[nsat] = -naxis2 ;
  sz = szz = 0.0 ;

  for (j=0; j<naxis2; j++) {
    for (i=0; i<naxis1; i++) {
      z = (int)(*(indata+i+j*naxis1)) ;
      if (z>cmax) {  // pixel is "saturated", tag it for obliteration
	jx = abs(isat[nsat] - i) ;
	jy = abs(jsat[nsat] - j) ;
	if ((jx>20) || (jy>20)) {
	  nsat++ ;
	  if (nsat > MAXSAT) {
	    printf("getsky: error, # saturated pixels > %d\n",MAXSAT);
	    exit(0);
	  }
	  isat[nsat] = i ;
	  jsat[nsat] = j ;
	}
      }
      if ((z<top) && (z>bot)) {  // pixel is in range, count it
	kk = (int)(aa * z + bb) ;
	*(hist+kk) += 1 ;
      }
      *(noise+i+j*naxis1) = z ; // fill array for FWHM calculation 
      sz += z ;
      szz += z*z ;
      npix++ ;
    }
  }

  // obliterate region around saturated stars.
  // Beware not to overwrite the array boundaries! [rwp/osu]

  if (nsat>0) {
    for (k=1; k<=nsat; k++) {
      j1 = jsat[k] - 20 ;
      j2 = jsat[k] + 20 ;
      i1 = isat[k] - 20 ;
      i2 = isat[k] + 20 ;
      for (j=j1; j<j2; j++) {
	for (i=i1; i<i2; i++) {
	  if (i < naxis1 && j < naxis2 && i > 0 && j > 0) {
	    *(noise+i+j*naxis1) = magic ;
	  }
	}
      }
    }
  }

  // find mode of histogram; set sky = mode 

  mx = 0 ;
  for (i=0; i<nhist; i++)
    if (*(hist+i) > mx) {
      jx = i ;
      mx = *(hist+i) ;
    }

  *sky = 1 + ((jx - bb) / aa) ;
  free(hist) ;
}

//---------------------------------------------------------------------------

#include <math.h>

#define max(a,b) ((a)>(b) ? (a) : (b))
          
void  
getfwhm(int naxis1,int naxis2,int ibor,float gain,int *noise,float *fwhm,
	float *tsmin, float *tsmax, float sky) 
{

  int numstr,j1,j2,i1,i2,i0,j0,z,zinmax,zmx ;
  int id,jd,i,j,zm,i1m,j1m,i2m,j2m,*nn,mm,ns ;
  float *sgx,*sgy,sg2x,sg2y,sxx,syy,fwhmx,fwhmy,sx,sy ;

  // use part of image (i1,j1) -> (i2,j2) 

  j1 = i1 = ibor - 1 ;
  j2 = naxis2 - ibor - 1 ;
  i2 = naxis1 - ibor - 1 ;
  numstr = 0 ;

  // allocate space for FWHM estimates for MAXSTARS stars 

  sgx = (float *)calloc(MAXSTARS,sizeof(float)) ;
  sgy = (float *)calloc(MAXSTARS,sizeof(float)) ;
  sg2x = sg2y = sxx = syy = 0.0 ;
  zinmax = -1 ;

  do {

    // find x,y with next highest count value - possible star, up to 20
    // total

    zmx = -1000000 ;
    for (j=j1; j<=j2; j++)
      for (i=i1; i<=i2; i++)
        if ((mm = *(noise+i+j*naxis1))<magic) {
          z = mm ;
          if (z>starmax) starmax = z ;
          if (z>zmx) {
            zmx = z ;
            i0 = i ;
            j0 = j ;
          }
        }
    zmx -= (int)(sky) ;
    if (zinmax<0) zinmax = zmx/10 ;
    if (zinmax<pstar*sky) zinmax = (int)(pstar*sky) ;

    zm = zmx / 2 ;
    
    // locate half-signal points in x and y for possible star centred at
    // i0,j0

    for (i=i0; i>i1; i--)
      if (*(noise+i+j0*naxis1)-sky < zm) break ;
    i1m = i ;

    for (i=i0; i<i2; i++)
      if (*(noise+i+j0*naxis1)-sky < zm) break ;
    i2m = i ;

    for (j=j0; j>j1; j--)
      if (*(noise+i0+j*naxis1)-sky < zm) break ;
    j1m = j ;

    for (j=j0; j<j2; j++)
      if (*(noise+i0+j*naxis1)-sky < zm) break ;
    j2m = j ;

    id = max(abs(i2m-i0),abs(i0-i1m)) ;
    jd = max(abs(j2m-j0),abs(j0-j1m)) ;
    id = max(id,jd) ;

    // id is estimate of full width 

    if (id>20) id = 20 ;
    nn = noise + (j0-id)*naxis1 ;

    // fit 2D Gaussian; failure -> fwhmx=0.0 

    fit2(i0,j0,id,naxis1,&fwhmy,&fwhmx,nn) ;

    // obliterate star to avoid finding it again 

    for (j=j0-4*id; j<j0+4*id; j++)
      for (i=i0-4*id; i<i0+4*id; i++)
	if (i < naxis1 && j < naxis2)
	  *(noise+i+j*naxis1) = magic;

    // store results for this object if valid and assemble stats 

    //    if (fwhmx > 0.01) {
    if (fwhmx > 0.01 && fwhmy > 0.01 && fwhmx < 15.0 && fwhmy < 15.0) {
      numstr++ ;
      *(sgx+numstr) = fwhmx ;
      *(sgy+numstr) = fwhmy ;
      sg2x +=fwhmx ;
      sg2y +=fwhmy ;
      sxx += fwhmx*fwhmx ;
      syy += fwhmy*fwhmy ;
    }

  } while((numstr<MAXSTARS) && (zmx>zinmax)) ;

  // NO STARS FOUND; test for this on return, but compute tsmin/tsmax 

  if (numstr<1) {
    printf("No suitable stars found on the image\n") ;
    *fwhm = -1.0 ;
    *tsmax = 0.8 *(10.0*zinmax - sky) ;
    if (*tsmax>starmax) *tsmax = 0.8 * starmax ;
    *tsmin = 5.0 * sqrt(sky/gain) ;
    return ;
  }

  // estimate upper and lower threshold for DoPhot 

  *tsmax = 0.8 *(10.0*zinmax - sky) ;
  if (*tsmax>starmax) *tsmax = 0.8 * starmax ;
  *tsmin = 5.0 * sqrt(sky/gain) ;

  // find mean, sigma for fwhm in both x and y 

  sg2x /= numstr ;
  sg2y /= numstr ;
  sxx /= numstr ;
  syy /= numstr ;
  sx = 2.0 * sqrt(sxx - sg2x*sg2x) ;
  sy = 2.0 * sqrt(syy - sg2y*sg2y) ;
  sxx = syy = 0.0 ;
  ns = 0 ;

  // reject outliers 

  for (i=1; i<numstr; i++) {
    j = 1 ;
    if (numstr>3)
      if (abs((int)(*(sgx+i)-sg2x))<sx && abs((int)(*(sgy+i)-sg2y))<sy)       
        j = 1 ;
      else
        j = 0 ;
    if (j==1) {
      sxx += *(sgx+i) ;
      syy += *(sgy+i) ;
      ns++ ;
    }
  }
  
  // recompute mean values 

  if (ns > 0) {
    fwhmx = sxx/ns ;
    fwhmy = syy/ns ;
    *fwhm = 0.5 * (fwhmx + fwhmy) ;
  }
  else {
    *fwhm = -1.0;
  }

}
