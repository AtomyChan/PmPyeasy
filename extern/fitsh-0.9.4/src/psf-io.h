/*****************************************************************************/
/* psf-io.h								     */
/*****************************************************************************/

#ifndef	__PSF_IO_H_INCLUDED
#define	__PSF_IO_H_INCLUDED	1

#include <stdio.h>		/* for type (FILE *) 			*/
#include <fits/fits.h>		/* for type (fits *) and (fitsimage *)	*/

#include "psf.h"		/* for type (psf *)			*/

/* psf_write():
   Dumps PSF 'p' in a plottable format into the file 'fw'. Plottable means
   can be plotted with `gnuplot`, with `splot 'psf.dat' u 1:2:3` or like...  */
int	psf_write(FILE *fw,psf *p);

/* psf_write_fits():
   Writes PSF 'p' into the file 'fw' in FITS format, as a NAXIS=3 image. 
   The third ('z') axis describes the polynomial variation of the PSF blocks.
   some extra headers (PSFSGRID, PSFHSIZE, PSFORDER, PSFOFFSX, PSFOFFSY, 
   PSFSCALE) are also written into the FITS image.			     */
int	psf_write_fits(FILE *fw,psf *p);

/* psf_parse_fits():
   Converts the FITS image 'img' into a PSF (stored in 'p'). The headers
   written by psf_write_fits() are assumed to exist and have consistent
   values (e.g. NAXIS1 should be equal to PSFSGRID*(2*PSFHSIZE+1), and 
   like that, see source and the related parts of the documentation.	     */
int	psf_parse_fits(fits *img,psf *p);

/* psf_read_fits():
   Reads the PSF (which is stored after in 'p') from the stream 'fr'. 
   This function uses psf_parse_fits() and the file should be a consistent
   PSF file (see constraints above).					     */
int	psf_read_fits(FILE *fr,psf *p);

#endif

