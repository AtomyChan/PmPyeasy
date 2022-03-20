//---------------------------------------------------------------------------
//
// fit2 - set up for 2D Gaussian fit to star image
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>

// Function Prototypes

void fit2(int ,int ,int ,int ,float *,float *,int *);
void free_matrix(float **,int ,int ,int ,int );
float **matrix(int ,int ,int ,int );
void bgauss(float **,int ,int ,int ,int ,int ,float ,float ,
	    int ,float *,float *);
void setsup(int mpt,float radg,float cen,float **f,float **sf,
	    float ra[4],int radflag);
float prodcr(float **a,int ja,float **b,int jb,int ns,int ne);
float prodrr(float **a,int ja,float **b,int jb,int ns,int ne);
void solve(float **a,float b[],int nx);
int gaures(float err[],float *xcen,float *ycen,float *sigx,float *sigy,
	   int xlen,int ylen);
void nrerror(char error_text[]);
float *vector(int nl,int nh);
int *ivector(int nl,int nh);
void free_vector(float *v,int nl,int nh);
void free_ivector(int *v,int nl,int nh);
void covsrt(float **covar,int ma,int lista[],int mfit);

// globals

float amag,ht,base,rms,radx,rady,texpos ;
int iter ;

void
fit2(int xc0,int yc0,int radius,int ncol,float *fwhmx,float *fwhmy,
     int *pdata)
{
  int *pdat, *p ;
  //  unsigned int num1,radflag,itlim ;
  int num1,radflag,itlim ;
  int plast,ninval,ix,iy ;
  float **dat,xc,yc ;
  float xx,yy ;
  char inout[30] ;

  num1 = (2*radius+1) ;
  dat = matrix(1,num1,1,num1) ;

  pdat = pdata + xc0 - radius - 1 ;
  plast = *pdat ;
  radflag = 0 ;  // find Gaussian radii. =1 to fit with fixed radii
  itlim = 10 ;
  ninval = 0 ;
  for(iy=1; iy<=num1; iy++) {
    p = pdat + (iy-1)*ncol  ;
    for(ix=1; ix<=num1; ix++) {
      p++ ;
      if(*p<65535) {
	dat[iy][ix] = (float) *p ;
	plast = *p ;
      }
      else {
	dat[iy][ix] = (float) plast ;
	ninval++ ;
      }
    }
  }
  bgauss(dat,num1,num1,num1,num1,radflag,xx,yy,itlim,&xc,&yc) ;
  *fwhmx = 1.665 * radx ;
  *fwhmy = 1.665 * rady ;
}

//---------------------------------------------------------------------------
// 
// fit_gauss() - 2-D least squares Gaussian fit
//
// routines:
//  bgauss - controls the fitting                         
//  solve  - solve symmetric set of simultaneous equations
//  setsup - set up terms for simultaneous equations      
//  prodrr - multiply row by row                          
//  prodcr - multiply column by row                       
//  gaures - updates fit and checks for convergence       
//

void
bgauss(float **apt,int xlen,int ylen,int mm,int nn,int radflag,
       float orix,float oriy,int itlim,float *xcen,float *ycen)
{
  float s[8],**d,e[7],**f,**g,**bpt ;
  float ra[4],rb[4],**sf,**sg ;
  static int ise[] = {0,6,1,2,4,3,5} ;
  static int isf[] = {0,4,1,2,1,3,1} ;
  static int isg[] = {0,4,1,1,2,1,3} ;
  float sum;
  int ne,nf,ng,itera,ninval,j,k,l ;
  double xmag, ymag ;
  char inout[20] ;

  f = matrix(1,xlen,1,4) ;
  g = matrix(1,xlen,1,4) ;
  d = matrix(1,6,1,6) ;
  sf = matrix(1,4,1,4) ;
  sg = matrix(1,4,1,4) ;
  bpt = apt ;
  *xcen = (float)mm / 2.0 + 0.5 ;
  *ycen = (float)nn / 2.0 + 0.5 ;
  if(radflag==0) {
    ne = 6 ;
    nf = 3 ;
    ng = 3 ;
    radx = mm/6.0 ;
    rady = nn/6.0 ;
  }
  else {
    ne = 4 ;
    nf = 2 ;
    ng = 2 ;
    radx = orix ;
    rady = oriy ;
  }
  amag = 0.0 ;
  ht = 0.001 ;
  base = 0.001 ;
  rms = 0.0 ;
  ninval = 0 ;
  sum = 0.0 ;
  for(k=1; k<=nn; k++)
    for(j=1; j<=mm; j++)
      sum += apt[k][j] ;
  
  itera = 0 ;
  iter = 0 ;
  do {
    setsup(mm,radx,*xcen,f,sf,ra,nf) ;
    setsup(nn,rady,*ycen,g,sg,rb,ng) ;
    s[6] = sum ;
    for(j=1; j<=mm; j++)
      f[j][4] = prodcr(apt,j,g,1,1,nn) ;
    for(j=1; j<=nn; j++)
      g[j][4] = prodrr(apt,j,f,1,1,mm) ;
    for(l=1; l<=nf; l++)
      s[l] = ra[l]*prodrr(f,4,f,l,1,mm) ;
    for(l=2; l<=ng; l++)
      s[l+2] = rb[l]*prodrr(g,4,g,l,1,nn) ;
    for(k=1; k<=ne; k++) {
      for(j=k; j<=ne; j++) {
	d[j][k] = sf[isf[k]][isf[j]]*sg[isg[k]][isg[j]] ;
	d[k][j] = d[j][k] ;
      }
      e[k] = s[ise[k]] ;
    }
    solve(d,e,ne) ;
    itera = gaures(e,xcen,ycen,&radx,&rady,mm,nn) ;
    iter += 1 ;
  } while(iter<itlim && itera==0) ;
  rms = 0.0 ;
  for(k=1; k<=nn; k++)
    for(j=1; j<=mm; j++) {
      sum = bpt[j][k] - ((ht)*f[j][1]*g[k][1] + (base)) ;
      rms += sum*sum ;
    }
  rms = (float) sqrt((double)rms/(nn*mm-1.0)) ;
  if(texpos<0.001)
    texpos = 1.0 ;
  xmag = (double)(ht) * (double)(radx) * (double)(rady) /(double)texpos ;
  if(xmag>0.0)
    ymag = log10(xmag) ;
  else
    ymag = 99.99 ;
  if(xmag > 0.0)                      //for approx agreement with fg_meas 
    amag = 18.75 - 2.5*(float)log10(xmag) ;// 20 - 2.5log10(pi) 
  else
    amag = 99.0 ;
  if(ht<=1.0e-5 || radx>50.0 || rady>50.0 || radx<1.2 || rady<1.2 ) {
    if(radflag==0) {
      radx = 0.0 ;
      rady = 0.0 ;
    }
    base = 0.0 ;
    ht = 0.0 ;
    rms = 0.0 ;
    amag = 50.0 ;
    *xcen = mm/2.0 + 0.5 ;
    *ycen = nn/2.0 + 0.5 ;
  }
  free_matrix(f,1,xlen,1,4) ;
  free_matrix(g,1,xlen,1,4) ;
  free_matrix(d,1,6,1,6) ;
  free_matrix(sf,1,4,1,4) ;
  free_matrix(sg,1,4,1,4) ;
}

#define sign(a) ((a)>=0 ? 1 : -1)

void
setsup(int mpt,float radg,float cen,float **f,float **sf,
       float ra[4],int radflag)
     //   int mpt,radflag ;
     //   float radg,cen,**f,**sf,ra[4] ;
{
  float rr,ss,zz,ww ;
  int j,k,l ;

  rr = 1.0/(radg) ;
  ra[1] = 1.0 ;
  ra[2] = 2.0*rr*rr ;
  ra[3] = rr*ra[2] ;
  for(j=1; j<=mpt; j++) {
    zz = j - (cen) ;
    ss = rr * zz ;
    ww = 0.0 ;
    if(fabs(ss)<6.0)
      ww = exp(-(double)ss*ss) ;
    f[j][1] = ww ;
    f[j][2] = zz * f[j][1] ;
    f[j][3] = zz * f[j][2] ;
  }
  for(k=1; k<=radflag; k++) {
    for(l=k; l<=radflag; l++) {
      sf[k][l] = ra[k]*ra[l]*prodrr(f,k,f,l,1,mpt) ;
      sf[l][k] = sf[k][l] ;
    }
    ww = 0.0 ;
    for(j=1; j<mpt; j++)
      ww += f[j][k] ;
    sf[4][k] = ww*ra[k] ;
    sf[k][4] = sf[4][k] ;
  }
  sf[4][4] = mpt ;
}

float 
prodcr(float **a,int ja,float **b,int jb,int ns,int ne)
{
  float sum ;
  int k ;
  
  sum = 0.0 ;
  for(k=ns; k<=ne; k++)
    sum += a[ja][k] * b[k][jb] ;
  return(sum) ;
}

float 
prodrr(float **a,int ja,float **b,int jb,int ns,int ne)
{
  float sum ;
  int k ;

  sum = 0.0 ;
  for(k=ns; k<=ne; k++)
    sum += a[k][ja] * b[k][jb] ;
  return(sum) ;
}

void
solve(float **a,float b[],int nx)
{
  int j,k ;
  float sum ;
  
  for (k=1; k<=nx; k++) {
    if (k!=1) {
      for(j=k; j<=nx; j++)
	a[j][k] -= prodcr(a,j,a,k,1,k-1) ;
    }
    if(k!=nx) {
      if(fabs(a[k][k])<1.0e-20)
	a[k][k] = sign(a[k][k]) * 1.0e-20 ;
      for(j=k+1; j<=nx; j++) {
	a[k][j] = a[j][k] ;
	a[j][k] = a[j][k]/a[k][k] ;
      }
    }
  }
  for(k=2; k<=nx; k++) {
    sum = 0.0 ;
    for(j=1; j<k; j++)
      sum += a[k][j] * b[j] ;
    b[k] -= sum ;
  }
  if(fabs(a[nx][nx])<1.0e-20)
    a[nx][nx] = sign(a[nx][nx]) * 1.0e-20 ;
  b[nx] /= a[nx][nx] ;
  if (nx!=1) {
    if (fabs(a[nx-1][nx-1])<1.0e-20)
      a[nx-1][nx-1] = sign(a[nx-1][nx-1]) * 1.0e-20 ;
    b[nx-1] = (b[nx-1] - a[nx-1][nx]*b[nx])/a[nx-1][nx-1] ;
    if(nx>=2)
      for(j=2; j<nx; j++) {
	sum = 0.0 ;
	for(k=1; k<=j; k++)
	  sum += a[nx-j][nx-j+k]*b[nx-j+k] ;
	if(fabs(a[nx-j][nx-j])<1.0e-20)
	  a[nx-j][nx-j] = sign(a[nx-j][nx-j])*1.0e-20 ;
	b[nx-j] = (b[nx-j]-sum)/a[nx-j][nx-j] ;
      }
  }
}

int
gaures(float err[],float *xcen,float *ycen,float *sigx,float *sigy,
       int xlen,int ylen)
{
  float rx0,ry0,hrx,hry,dx,dy,drx,dry ;

  rx0 = fabs(*sigx) ;
  ry0 = fabs(*sigy) ;
  hrx = rx0/2.0 ;
  hry = ry0/2.0 ;
  ht = err[2] ;
  if (fabs(ht)<1.0e-6)
    ht = sign(ht) * 1.0e-6 ;
  dx = err[3]/(ht) ;
  dy = err[4]/(ht) ;
  if (ht>0.0) {
    drx = err[5]/(ht) ;
    dry = err[6]/(ht) ;
  }
  else {
    dx *= -1.0 ;
    dy *= -1.0 ;
    drx = 0.0 ;
    dry = 0.0 ;
  }
  base = err[1] ;
  if(fabs(dx)>hrx)
    dx = sign(hrx) * dx ;
  if(fabs(dy)>hry)
    dy = sign(hry) * dy ;
  *xcen += dx ;
  if(*xcen<3.0)
    *xcen = 3.0 ;
  else if (*xcen > xlen-2.0)
    *xcen = xlen - 2.0 ;
  *ycen += dy ;
  if(*ycen<3.0)
    *ycen = 3.0 ;
  else if (*ycen > ylen-2.0)
    *ycen = ylen - 2.0 ;
  *sigx = rx0 + drx ;
  *sigy = ry0 + dry ;
  if(*sigx<0.5)
    *sigx = 0.5 ;
  if(*sigy<0.5)
    *sigy = 0.5 ;
  if(fabs(dx)>0.1 || fabs(dy)>0.1 || fabs(drx)>0.005 || fabs(dry)>0.005)
    return(0) ;
  else
    return(1) ;
}

//----------------------------------------------------------------
//
//  Vector and matrix stuff       
//

void 
nrerror(char error_text[])
{

  fprintf(stderr,"Runtime error...\n") ;
  fprintf(stderr,"%s\n",error_text) ;
  fprintf(stderr,"...now exiting to system...\n") ;
  fprintf(stderr,"(program: fwhmsky)\n");
  exit(1) ;
}

float 
*vector(int nl,int nh)
{
  float *v ;

  v = (float *)malloc((unsigned)(nh-nl+1)*sizeof(float)) ;
  if(!v)
    nrerror("Allocation failure in vector()") ;
  return(v-nl) ;
}

int 
*ivector(int nl,int nh)
{
  int *v ;

  v = (int *)malloc((unsigned)(nh-nl+1)*sizeof(int)) ;
  if(!v)
    nrerror("Allocation failure in ivector()") ;
  return(v-nl) ;
}

void 
free_vector(float *v,int nl,int nh)
{
  free((char*)(v+nl)) ;
}


void 
free_ivector(int *v,int nl,int nh)
{
  free((char*)(v+nl)) ;
}

void 
covsrt(float **covar,int ma,int lista[],int mfit)
{
  int i,j ;
  float swap ;

  for(j=1; j<ma; j++)
    for(i=j+1; i<=ma; i++)
      covar[i][j] = 0.0 ;
  for(i=1; i<mfit; i++)
    for(j=i+1; j<=mfit; j++) {
      if(lista[j] > lista[i])
	covar[lista[j]][lista[i]] = covar[i][j] ;
      else
	covar[lista[i]][lista[j]] = covar[i][j] ;
    }
  swap = covar[1][1] ;
  for(j=1; j<=ma; j++) {
    covar[1][j] = covar[j][j] ;
    covar[j][j] = 0.0 ;
  }
  covar[lista[1]][lista[1]] = swap ;
  for(j=2; j<=mfit; j++)
    covar[lista[j]][lista[j]] = covar[1][j] ;
  for(j=2; j<=ma; j++)
    for(i=1; i<=j-1; i++)
      covar[i][j] = covar[j][i] ;
}


float 
**matrix(int nrl,int nrh,int ncl,int nch)
{
  int i ;
  float **m ;

  m = (float **)malloc((unsigned)(nrh-nrl+1)*sizeof(float*)) ;
  if(!m)
    nrerror("allocation failure 1 in matrix()") ;
  m -=nrl ;
  for(i=nrl; i<=nrh; i++) {
    m[i] = (float *)malloc((unsigned)(nch-ncl+1)*sizeof(float)) ;
    if(!m[i])
      nrerror("allocation failure 2 in matrix()") ;
    m[i] -= ncl ;
  }
  return(m) ;
}

void 
free_matrix(float **m,int nrl,int nrh,int ncl,int nch)
{
  int i ;

  for(i=nrh; i>=nrl; i--)
    free((char *)(m[i]+ncl)) ;
  free((char *)(m+nrl)) ;
}
