/*****************************************************************************/
/* psf-base.c								     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Functions related to basic/common PSF handling. Actually, there is only   */
/* one function here: psf_symmetrize()...		 	             */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* (c) 2006; Pal, A. (apal@szofi.elte.hu)			             */
/*****************************************************************************/

#include "psf.h"
#include "psf-base.h"

/*****************************************************************************/

int psf_symmetrize(psf *pdata)
{
 double	**parr,p,***psfstack;
 int	i,j,w,nrd,nval;

 w=pdata->grid*(pdata->hsize*2+1);
 nval=(pdata->order+1)*(pdata->order+2)/2;

 psfstack=pdata->coeff;

 for ( nrd=nval ; nrd>0 ; nrd--,psfstack++ )
  {	parr=*psfstack;
	for ( i=0 ; i<(w+1)/2 ; i++ )
	 {	for ( j=i ; j<(w+1)/2 ; j++ )
		 {	p=parr[i][j]+parr[j][i]+
			  parr[w-1-i][j]+parr[w-1-j][i]+
			  parr[i][w-1-j]+parr[j][w-1-i]+
			  parr[w-1-i][w-1-j]+parr[w-1-j][w-1-i];
			p=p/8.0;
			parr[i][j]=parr[j][i]=
			parr[w-1-i][j]=parr[w-1-j][i]=
			parr[i][w-1-j]=parr[j][w-1-i]=
			parr[w-1-i][w-1-j]=parr[w-1-j][w-1-i]=p;
		 }
	 }
  }
 return(0);	
}

/*****************************************************************************/

