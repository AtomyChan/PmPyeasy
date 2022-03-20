/*****************************************************************************/
/* fits-common.h							     */
/*****************************************************************************/

#ifndef	__FITS_COMMON_H_INCLUDED
#define	__FITS_COMMON_H_INCLUDED	1

#ifndef	_FITS_SOURCE
#error  "This file (fits-common.h) is for internal usage, do not include directly."
#endif

#include <fits/fits.h>

#ifdef		FITS_TAPE_BLOCKSIZE
#undef		FITS_TAPE_BLOCKSIZE
#endif
#define		FITS_TAPE_BLOCKSIZE		2880

#ifdef		FITS_TAPE_CARDIMAGESIZE
#undef		FITS_TAPE_CARDIMAGESIZE
#endif
#define		FITS_TAPE_CARDIMAGESIZE		  80

#define		FITS_TAPE_CARDIMAGECOUNT	  36	

/*****************************************************************************/

#define		BUFSIZE		16384
#define		HDRBLOCKS	64

/*****************************************************************************/

#define		HDR_NAXIS	0
#define		HDR_NAXIS1	1
#define		HDR_NAXIS2	2
#define		HDR_NAXIS3	3
#define		HDR_NAXIS4	4
#define		HDR_BITPIX	5
#define		HDR_BSCALE	6
#define		HDR_BZERO	7
#define		HDR_EXTEND	8
#define		HDR_XTENSION	9
#define		HDR_INHERIT	10
#define		HDR_SIMPLE	11
#define		HDR_ORIGIN	12
#define		HDR_GAIN	13
#define		HDR_END		14

#define		HDR_TFIELDS	15
#define		HDR_TBCOL	16
#define		HDR_TFORM	17
#define		HDR_TSCAL	18
#define		HDR_TZERO	19
#define		HDR_TNULL	20
#define		HDR_TTYPE	21
#define		HDR_TUNIT	22
#define		HDR_PCOUNT	23
#define		HDR_GCOUNT	24

extern	char	*headers[];
extern	char	*comment_fits_standard;
extern	int	substantial_headers[];

/*****************************************************************************/

void *	fits_tensor_alloc_arr(int typesize,int rank,int *arr);
int	fits_tensor_free(void *tensor);

/*****************************************************************************/

int	fits_arch_is_swapped(void);
int	fits_swap_line_bytes(unsigned char *wbuff,int bs,int sx);

/*****************************************************************************/

int	fits_cb_read  (void *param,void *dst,int length);
int	fits_cb_write (void *param,void *src,int length);
int	fits_cb_skip  (void *param,int length);
int	fits_cb_is_end(void *param);

typedef struct
 {	char	*buffer;
	int	length;
 } fitsmemwrite;

int	fits_cb_mem_write(void *param,void *src,int length);

/*****************************************************************************/

#endif
                
/*****************************************************************************/
                                         
                 
