/*****************************************************************************/
/* fitsmask.h								     */
/*****************************************************************************/

#ifndef	__FITS_XTN_H_INCLUDED
#define	__FITS_XTN_H_INCLUDED	1

/*****************************************************************************/

#include <fits/fits.h>

/*****************************************************************************/

/* fits_mask_read_from_header():
   Creates a new mask using the header information (see 'MASKINFO' headers
   or any other optional headers can be set by 'hdr' if it is not NULL)
   of the FITS image 'img'. If 'sx' and 'sy' are valid sizes (thus,
   both of them are positive), the mask returned will have a size of
   'sx' times 'sy', otherwise the size of the mask will be set using
   the appropriate fields of 'img'. If the creation of the mask fails
   (due to any reason, allocation error, or these size fields are invalid)
   the function returns NULL, otherwise it returns the pointer of the
   newly allocated mask (which can be free with fits_mask_free()).	     */
char **	fits_mask_read_from_header(fitsheaderset *header,int sx,int sy,char *hdr);

/* fits_mask_mask_from_header():
   Performs a logical 'and' operation between the mask 'mask' and the
   mask of the image 'img' (if there is any mask described by 'mh' or
   MASKINFO headers). The result is also stored in 'mask' (so the 
   original mask 'mask' is going to be overwritten by the new one).
   If there's no mask information in 'img', the mask 'mask' is unchanged.    */
int fits_mask_mask_from_header(char **mask,fitsheaderset *header,int sx,int sy,char *mh);

int fits_mask_mask_line(char *line,fitsheaderset *header,int sx,int sy,int y0,char *maskhdr);
int fits_mask_mask_more_line(char **lines,fitsheaderset *header,int sx,int sy,int y0,int ny,char *maskhdr);

/* fits_mask_duplicate():
   Duplicates the mask 'mask' (which has a size of 'sx' by 'sy'). If
   the duplication was not successful (usually due to insufficient memory),
   the function returns NULL, otherwise it returns the new mask.	     */
char ** fits_mask_duplicate(char **mask,int sx,int sy);

/* fits_mask_export_as_header():
   Compresses the mask 'mask' (which has a size of 'sx' times 'sy' or
   if these values are nonpositive, the size is derived using 'img')
   and stores it in the FITS header area of 'img'. The storage of the mask
   is almos human-readable (ascii numbers are stored), it's not so 
   important to talk about the algorithm (which is quite primitive, see
   the source code for more). If the flag 'cl' is true, the function clears
   all mask information from the image 'img' before store the new mask.
   It's highly recommended to set 'cl' to true. If the function fails,
   it returns a nonzero value, otherwise it returns 0. If 'hdr' is not NULL,
   the header name specified by it is used, otherwise the default 'MASKINFO'
   header is set by the function.					     */
int	fits_mask_export_as_header(fitsheaderset *header,int cl,char **mask,int sx,int sy,char *hdr);

/* fits_mask_create_empty():
   Creates an empty mask, with the size of 'sx' by 'sy' (if they 
   are positive, otherwise the size if derived from 'img' if it is not NULL).
   The mask has all true values (actually filled with 1's). If the
   creation fails (due to any reason, allocation error, the 'sx'/'sy' values 
   are invalid) the function returns NULL, otherwise it returns the pointer of 
   the newly allocated mask (which can be freed with fits_mask_free()).      */
char ** fits_mask_create_empty(int sx,int sy);

/* fits_mask_expand_false():
   Expands the false area in the mask 'mask' with a size of 'sx'
   by 'sy' with a border with a size of 'dis'. Note that the outbound area
   of the matrix 'cmin' is treated as a false area. The newly created matrix 
   is returned. If NULL is returned, the allocation of the new matrix failed.
   The primarily usage of the function is to determine where can a point
   of an image be convolved with a given kernel. If the matrix 'cmin' contains
   the valid pixel mask of an image (with a size of 'sx' by 'sy') and 'dis' is 
   the half-size of the kernel used for the convolution, the returned array
   contains the mask which has true elements where the convolution can be 
   performed, so the kernel used doesn't intersect into any false pixel or into
   the outbound of the image.						     */
char **	fits_mask_expand_false(char **mask,int sx,int sy,int dis,int imv,int smv,int expand_border);

int	fits_mask_expand_logic(char **mask,int sx,int sy,int dis,int imv,int smv,int expand_border);

/* fits_mask_and(), fits_mask_or():
   Performs a logical 'and' or 'or' operation between the mask 'mask' and
   'other' (both of them has a size of 'sx' by 'sy'). The result is 
   stored in 'mask' (the original mask 'mask' is going to be overwritten).
   If 'other' is NULL (or its all elements are true), the mask 'mask'
   is unchanged.							     */
int	fits_mask_and(char **mask,int sx,int sy,char **other);
int	fits_mask_or(char **mask,int sx,int sy,char **other);

/* fits_mask_create_floyd():
   Creates a Floyd-Steinberg dithered mask with the fill ratio of 'a'/'b'.
   The other mask creation parameters are like to fits_mask_create_empty().  */
char ** fits_mask_create_floyd(int sx,int sy,int a,int b,int smv);

/* fits_mask_mask_nans():
   Marks the pixels in the mask 'mask' where img->data is not a number 
   (truly, ! isfinite()) with the flag 'setmask'.			     */
int	fits_mask_mark_nans(fitsimage *img,char **mask,int setmask);

/* fits_mask_trim():
   Trims the mask 'mask' (with the original size of 'sx' times 'sy'). 	     */
int	fits_mask_trim(char ***mask,int sx,int sy,int x0,int y0,int nx,int ny,int outer);

/* fits_mask_is_clean():
   This function returns true (nonzero value) if all of the mask elements
   are zero, i.e. the mask/image is practically 'clean' from nasty pixels.   */
int	fits_mask_is_clean(char **mask,int sx,int sy);

/* fits_mask_free():
   Releases the memory are of 'mask' if it is not NULL. The function
   always return 0.							     */
int	fits_mask_free(char **mask);

/*****************************************************************************/

#endif

/*****************************************************************************/
                                                             
