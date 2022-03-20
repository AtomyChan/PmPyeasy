/*****************************************************************************/
/* fits.h 								     */
/*****************************************************************************/

#ifndef	__FITS_H_INCLUDED
#define	__FITS_H_INCLUDED	1

/*****************************************************************************/

#include <stdio.h>			/* (FILE *)	*/

/*****************************************************************************/

/* These values are highly parts of the FITS standard, never change them */
#define		FITS_TAPE_BLOCKSIZE		2880
#define		FITS_TAPE_CARDIMAGESIZE		  80

/*****************************************************************************/

#define		FITS_MAX_NAXIS	17	/* maximum image dimension	 */
					/* if you need more, change this */
					/* and recompile _all_ stuff.	 */

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#define		FITS_EMPTY	0	/* Empty (80 spaces)		*/
#define		FITS_VSTR	1	/* Normal string, 'quoted'...	*/
#define		FITS_VBOOLEAN	2	/* Boolean (stored as int)	*/
					/* e.g: SIMPLE = T, EXTEND = F	*/
#define		FITS_VINT	3	/* Integer (stored as int)	*/
#define		FITS_VDOUBLE	4	/* Real (stored as double)	*/
#define		FITS_VCOMMENT	5	/* Commentary text w/wo keyword	*/

#define		FITS_SH_FIRST		0	/* See fits_set_header_any() */
#define		FITS_SH_LAST		1
#define		FITS_SH_ADD		2
#define		FITS_SH_INSERT		3
#define		FITS_SH_FORCEFIRST	4
#define		FITS_SH_BEGIN		5

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#define		FITSERR_OPEN		1
#define		FITSERR_ALLOC		2
#define		FITSERR_FRAME_MISSING	3
#define		FITSERR_FRAME_AMBIG	4
#define		FITSERR_FRAME_INVALID	5
#define		FITSERR_IMAGE		6
#define		FITSERR_SCALE		7

/******************************************************************************
The FITS header element sturcure: describes one header element:
 - name: the keyword of the card image (in the standard FITS files, it should
   not have a lenght more than 8 characters (excluding the trailing zero);
 - comment: the comment for the given image (the string after the first '/';
 - type: the type of the value stored in the given card image, see the 
   declarations of FITS_EMPTY, ...;
 - vstr: if type==FITS_VSTR or FITS_VCOMMENT, the string is stored in vstr
   (including the trailing zero and excluding the single quotation marks);
 - vint: if type==FITS_VINT or FITS_VBOOLEAN the data is here;
 - vdouble: if type==FITS_VDOUBLE, the data is here.
The FITS header structure: describes a whole header:
 - hdrs: array of the header elements
 - nhdr: number of the header elements (trailing 'END' is not included!)
 - ahdr: the allocated memory for headers (should not be less than 'nhdr')
******************************************************************************/

typedef struct
 {	char		name[FITS_TAPE_CARDIMAGESIZE];
	char		comment[FITS_TAPE_CARDIMAGESIZE];
	char		vstr[FITS_TAPE_CARDIMAGESIZE];
	int		vtype;
	int		vint;
	double		vdouble;
 } fitsheader;

typedef struct			/* FITS main/extension header information:   */
 {	fitsheader	*hdrs;	/*  - header array ('END' is not included!)  */
	int		nhdr;	/*  - number of headers 		     */
	int		ahdr;   /*  - currently allocated entries for 'hdrs' */
 } fitsheaderset;

/*****************************************************************************/

typedef struct
 {	double	bscale;
	double	bzero;
 } fitsscale;

typedef struct
 {	int		sx,sy;	/*	The fits file contanis an image:     */
	int		bit;	/*	with a size 'sx' x 'sy', readed at   */
	double		**data;	/*	the depth 'bit' and stored in 'data' */

	int		dim;			/* if NAXIS != 2...	     */
	int		naxis[FITS_MAX_NAXIS];	/* values of NAXIS[1-4]	     */
	void		*vdata;

	void		*allocdata;	/* alloc. image (freed if not NULL)  */

	fitsscale	curr,read;	/* scaling and bias values	     */
 } fitsimage;

typedef struct				/* Fields: ASCII (Text) table	     */
 {	int		colindex;	/* value of TBCOLn - 1 		     */
	char		format[12];	/* value of TFORMn		     */
	fitsscale	scale;		/* values of TSCALn and TZEROn	     */
	char		type[32];	/* value of TTYPEn		     */
	char		unit[32];	/* value of TUNITn		     */
	char		null[32];	/* value of TNULLn		     */
 } fitstfield;

typedef struct
 {	int		nrow;		/* value of NAXIS2		     */
	int		rowsize;	/* value of NAXIS1		     */
	int		ntfield;	/* value of TFIELDS		     */
	fitstfield	*tfields;
	unsigned char	**data;
	void		*allocdata;
 } fitsttable;

/* Definitions of valid TFORM data types in BINTABLE extensions */
/* See [1], p43., Table 8.5					*/
#define		FTF_LOGICAL		('L')
#define		FTF_BIT			('X')
#define		FTF_BYTE		('B')
#define		FTF_UINT8		FTF_BYTE
#define		FTF_SHORT		('I')
#define		FTF_INT16		FTF_SHORT
#define		FTF_LONG		('J')
#define		FTF_INT32		FTF_LONG
#define		FTF_LONGLONG		('K')
#define		FTF_INT64		FTF_LONGLONG
#define		FTF_STRING		('A')
#define		FTF_CHARACTER		FTF_STRING
#define		FTF_CHAR		FTF_STRING
#define		FTF_FLOAT		('E')
#define		FTF_REAL32		FTF_FLOAT
#define		FTF_DOUBLE		('D')
#define		FTF_REAL64		FTF_DOUBLE
#define		FTF_FLOATCOMPLEX	('C')
#define		FTF_COMPLEX32		FTF_FLOATCOMPLEX
#define		FTF_DOUBLECOMPLEX	('M')
#define		FTF_COMPLEX64		FTF_DOUBLECOMPLEX
#define		FTF_ARRAY		('P')

typedef struct				/* Fields: binary table		     */
 {	int		form;		/* derived from TFORMn (see FTF_*)   */
	int		repeat;		/* derived from TFORMn		     */
	int		basesize;	/* derived from 'form'		     */
	int		offset;		/* derived from previous values      */
	fitsscale	scale;		/* values of TSCALn and TZEROn	     */
	char		type[32];	/* value of TTYPEn		     */
	char		unit[32];	/* value of TUNITn		     */
	char		null[32];	/* value of TNULLn		     */
 } fitsbfield;

typedef struct
 {	int		nrow;		/* value of NAXIS2		     */
	int		rowsize;	/* value of NAXIS1		     */
	int		nbfield;	/* value of TFIELDS		     */
	fitsbfield	*bfields;
	int		pcount;		/* value of PCOUNT		     */
	unsigned char	**data;
	void		*allocdata;
 } fitsbtable;

typedef union
 {	fitsimage	i;
	fitsttable	t;
	fitsbtable	b;
 } fitsextensiondata;

#define		FITS_EXT_IMAGE		1
#define		FITS_EXT_TABLE		2
#define		FITS_EXT_BINTABLE	3

typedef struct
 {	int			type;
	fitsheaderset		header;
	fitsextensiondata	x;
 } fitsextension;

/******************************************************************************
The complete FITS data:
 - (fitsheaderset)header: the complete primary header of the FITS data
 - (fitsimage)i:	  the primary image, if any (so i.vdata != NULL).
 - (fitsextension[])xtns: array of extensions
 - (int)nxtn:		  number of extensions
 - length: if the FITS contains something else than an image (so, the 
   header doesn't have appropriate NAXIS, BITPIX, NAXIS1 and NAXIS2 fields)
 - rawdata: the data in the FITS file, a dynamically allocated array
   of 'length' bytes
******************************************************************************/

typedef struct
 {	fitsheaderset	header;
	fitsimage	i;	
	fitsextension	*xtns;
	int		nxtn;
	int		length;		
	char		*rawdata;	
 } fits;

/*****************************************************************************/
/* Function prototypes							     */
/*****************************************************************************/

/******************************************************************************
 Note that I/O routines with trailing `_cb` in their name has the same 
 functionality as the originating functions with the appropriate file 
 stream (FILE *fr or FILE *fw), but these functions - instead of reading from
 or write to file streams - are using callback functions for reading / writing.

 The function cb_read(void *param,void *dst,int length) is for reading 'length'
 bytes to the buffer 'dst', while the optional parameter 'param' is passed 
 directy at each call. This function should return a negative number on error, 
 0 when end-of-input-data is reached, otherwise the total number of bytes read.
 If 'dst' is NULL, the function _must_ jump over the next 'length' bytes of the
 input data and should not fail (therefore, in this case the callback should
 behave like lseek() and fseek() calls with SEEK_CUR argument; but here
 the value of 'length' is always non-negative).

 The function cb_write(void *param,void *src,int length) is for writing
 'length' bytes from the buffer 'src', while the optional parameter 'param'
 is passed directly at each call. This function should return a negative
 integer on error, otherwise the total number of bytes written. If 'src' is 
 NULL, the function can do anything, e.g. it can return 0 (0 bytes has been
 written from NULL). 

 For an impementation, see fits_cb_read() and fits_cb_write() in the module
 fits-common.c. These functions do reading/writing from file streams:
 the parameter 'param' is a (void *)-casted form of (FILE *) and these
 functions cal fread(), fseek() and fwrite() directly. 

 This whole callback-stuff is the abstraction of the underlying I/O of FITS
 data handling. In practice, it can also be used for reading/writing FITS data
 from/to network sockets, pipes and also for dumping the whole FITS data
 into a machine-independent format (especially in the native FITS binary 
 representation) in the conventional memory.
******************************************************************************/

/* fits_create():
   Creates an empty FITS descriptor (with no header, no data and no rawdata).*/
fits *		fits_create(void);
/* fits_duplicate():
   Duplicates the FITS 'img' (even when it is not an image).                 */
fits *		fits_duplicate(fits *img);
/* fits_duplicate_empty():
   Duplicates the FITS 'img' and if it was an image, fills the array 'data' 
   with zeros.								     */
fits *		fits_duplicate_empty(fits *img);
/* fits_free_image():
   Frees the array 'data' of the FITS 'img', if it was an image. After the
   release of the array, it sets the 'data' to NULL and resets 'sx' and 'sy'.*/
void		fits_free_image(fits *img);
/* fits_free():
   Frees the whole FITS 'img'.                                               */
void		fits_free(fits *img);

/* fits_tapeblock_size():
   Returns the required size to store 'size' bytes in FITS tape blocks.
   Practically, this function round the value of 'size' up to the nearest
   multiple of FITS_TAPE_BLOCKSIZE, which is 2880.			     */
int		fits_tapeblock_size(int size);

/* fits-header.c */ /*********************************************************/

/* fits_headerset_reset():
   Resets the headerset data (set all values to NULL and 0). Must be called
   before any fits_headerset_...() call if local automatic header structures
   want to be used. The function doesn't do any memory release.		     */
int		fits_headerset_reset(fitsheaderset *header);
/* fits_headerset_duplicate():
   Duplicates the headerset 'act'.					     */
int		fits_headerset_duplicate(fitsheaderset *r,fitsheaderset *act);
/* fits_headerset_free():
   Releases and resets the headerset data 'header'.			     */
int		fits_headerset_free(fitsheaderset *header);

/* fits_headerset_read(), fits_headerset_write(),
   fits_headerset_read_cb(), fits_headerset_write_cb():
   Reads/writes the complete headerset from/to 'fr'/'fw' or via callback.    */
int		fits_headerset_read(FILE *fr,fitsheaderset *header);
int		fits_headerset_write(FILE *fw,fitsheaderset *header);

int		fits_headerset_read_cb(int (*cb_read)(void *,void *,int),
			void *param,fitsheaderset *header);
int		fits_headerset_write_cb(int (*cb_write)(void *,void *,int),
			void *param,fitsheaderset *header);

/* fits_headerset_get_count():
   Counts how many header elements in the headerset have a name of 'hdr'.
   The comparison is case-sensitive, so, according to the FITS standards,
   the 'hdr' should be an uppercase-string.                                  */
int		fits_headerset_get_count(fitsheaderset *header,char *hdr);

/* fits_headerset_get_id():
   Returns the array index of the 'cnt'th copy of the header 'hdr' in
   the headerset 'header'. The 'cnt' should be between zero
   (the first such header) and fits_headerset_get_count(header,hdr)-1 (the 
   last such header). Otherwise, or if the header 'hdr' doesn't exist 
   in the headerset, the function returns a negative value.		     */
int		fits_headerset_get_id(fitsheaderset *header,char *hdr,int cnt);
/* fits_headerset_get_header(), fits_headerset_get_uniq_header():
   Returns the pointer to the header structure of 'cnt'th copy of the header 
   'hdr' in the headerset 'header'. This function does almost the same what 
   the fits_headerset_get_id does, but returns the pointer itself instead of
   the array index. If the header 'hdr' doesn't exist in the headerset or 'cnt'
   is negative or larger than fits_headerset_get_count(), it returns NULL. 
   The function fits_headerset_get_uniq_header() works like above but assumes
   the header element named 'hdr' to be unique. So if there is more than
   one header element with this name, this function also returns NULL.       */
fitsheader *	fits_headerset_get_header(fitsheaderset *header,
			char *hdr,int cnt);
fitsheader *	fits_headerset_get_uniq_header(fitsheaderset *header,
			char *hdr);


/* fits_headerset_get_as_double():
   Reads the numerical value from the first occurence of the header element
   named "hdr". If 'is_ambigous_allowed' is false (zero), just one copy of 
   the header 'hdr' is allowed (see also fits_headerset_get_uniq_header()). 
   Otherwise the function interprets the first occurence of the header element
   named "hdr". If the type of the header field with such name is an integer 
   or real (consequently, numerical), the function converts it to double, 
   stores it in  'ret' and returns zero. In all other cases -- there's no 
   such header element named "hdr", more than one "hdr" exist if 
   'is_ambigous_allowed' is false, or the header value is not numerical -- the 
   function returns a nonzero value.					     */
int		fits_headerset_get_as_double(fitsheaderset *header,
		char *hdr,double *ret,int is_ambigous_allowed);

/* fits_headerset_append(), fits_headerset_insert():
   Low-level functions to append or insert a single header element after or 
   before the headerset 'header'. The pointer to the newly created header 
   element is returned (or NULL if any error occured).			     */
fitsheader *	fits_headerset_append(fitsheaderset *header);
fitsheader *	fits_headerset_insert(fitsheaderset *header);

/* fits_headerset_set_any():
   Creates or updates a header element with a name of 'hdr' in the headerset.
   The behaviour of this function depends on the rule 'rule':
   - if 'rule' is FITS_SH_FIRST, the function returns a pointer to the first 
     element with the name of 'hdr' or creates a new header field with the name
     'hdr' if it didn't exist before the call (and then, the returned pointer 
     points to this newly created header field),
   - if 'rule' is FITS_SH_LAST, the function returns a pointer to the last 
     occurence of the header with the name of 'hdr' or creates a new header 
     field with the name  'hdr' if it didn't exist before the call,
   - if 'rule' is FITS_SH_ADD, the function creates a new header field with the
     name 'hdr' and the returned pointer points to this newly created header 
     field.
   In all cases, the comment field of the returned header is updated to 
   'c' if it was non-NULL (otherwise, the comment of the header is not
   changed or set to NULL if it is a new header). However, this function might 
   have had to be written as a static function (used by the other 
   fits_set_header_*() functions), sometimes it can be useful if low-level or 
   automatic manipulation of the headers is needed in any application. Use it 
   carefully.								     */
fitsheader *	fits_headerset_set_any(fitsheaderset *header,
			char *hdr,int rule,char *comment);

/* fits_headerset_set_[integer,double,string,boolean]():
   These four routines create or update a header element with the name of "hdr"
   in the headerset 'header'. The rules described by 'rule' are the same in 
   fits_headerset_set_any(). The structure element vtype is updated (to 
   FITS_VINT, FITS_VDOUBLE, FITS_VSTR or FITS_VBOOLEAN, respectively), and 
   the new value of the field is set. If 'c' is non-NULL, the comment of 
   the given header is updated also.					     */
int		fits_headerset_set_integer(fitsheaderset *header,
			char *hdr,int rule,int val   ,char *comment);
int		fits_headerset_set_double (fitsheaderset *header,
			char *hdr,int rule,double val,char *comment);
int		fits_headerset_set_string (fitsheaderset *header,
			char *hdr,int rule,char *str ,char *comment);
int		fits_headerset_set_boolean(fitsheaderset *header,
			char *hdr,int rule,int vbool ,char *comment);

/* fits_headerset_delete(), fits_headerset_delete_all():
   The function fits_headerset_delete() deletes the 'k'th occurrence of the 
   header element named "hdr" from the headerset 'header'. The function
   fits_headerset_delete_all() deletes all occurrences of header elements
   named "hdr" from the headerset 'header'.				     */
int		fits_headerset_delete(fitsheaderset *header,char *hdr,int k);
int		fits_headerset_delete_all(fitsheaderset *header,char *hdr);

/* fits_headerset_copy(), fits_headerset_merge():
   Copies or merges headersets (from 'h2' to 'h1' and from 'h' to 'xh').     */
int		fits_headerset_copy(fitsheaderset *h1,fitsheaderset *h2);
int		fits_headerset_merge(fitsheaderset *xh,fitsheaderset *h,
			int is_inherit);

/* fits_headerset_is_extension():
   Checks if the headerset 'header' is an extension header unit or not. If it
   is not an extension header (so, primary), the function returns 0, otherwise
   it returns a positive value indicating the type of the extension 
   (see the definitions: FITS_EXT_IMAGE, FITS_EXT_TABLE, FITS_EXT_BINTABLE). */
int		fits_headerset_is_extension(fitsheaderset *header);

/* fits-table.c */ /**********************************************************/

int		fits_table_get_params(fitsheaderset *header,fitsttable *ft);
int		fits_table_alloc(fitsttable *ft);
int		fits_table_read(FILE *fr,fitsttable *ft);
int		fits_table_read_cb(int (*cb_read)(void *,void *,int),
			void *param,fitsttable *ft);
int		fits_table_skip(FILE *fr,fitsttable *ft);
int		fits_table_skip_cb(int (*cb_read)(void *,void *,int),
			void *param,fitsttable *ft);

int		fits_table_set_params(fitsheaderset *header,fitsttable *ft);
int		fits_table_free(fitsttable *ft);

int		fits_table_write(FILE *fw,fitsttable *ft,int is_pad);
int		fits_table_write_cb(int (*cb_read)(void *,void *,int),
			void *param,fitsttable *ft,int is_pad);

int		fits_table_duplicate(fitsttable *ret,fitsttable *src,int flag);

/* fits-bintable.c */ /*******************************************************/

int		fits_bintable_form_basesize(int form);
char *		fits_bintable_form_cname(int form);
int		fits_bintable_get_params(fitsheaderset *header,fitsbtable *fb);

int		fits_bintable_create_fields(fitsbtable *fb,int,int nbf,...);
int		fits_bintable_check_fields(fitsbtable *fb,int nbf,...);
int		fits_bintable_set_xtr_params(fitsbtable *fb,int n,
			char *type,char *unit,char *null);
int		fits_bintable_set_xtr_scale(fitsbtable *fb,int n,
			double bscale,double bzero);

int		fits_bintable_alloc(fitsbtable *fb);
int		fits_bintable_read(FILE *fr,fitsbtable *fb);
int		fits_bintable_read_cb(int (*cb_read)(void *,void *,int),
			void *param,fitsbtable *fb);
int		fits_bintable_skip(FILE *fr,fitsbtable *fb);
int		fits_bintable_skip_cb(int (*cb_read)(void *,void *,int),
			void *param,fitsbtable *fb);

int		fits_bintable_set_offsets(fitsbtable *fb);
int		fits_bintable_set_params(fitsheaderset *header,fitsbtable *fb);

int		fits_bintable_free(fitsbtable *fb);

int		fits_bintable_write(FILE *fw,fitsbtable *fb,int is_pad);
int		fits_bintable_write_cb(int (*cb_read)(void *,void *,int),
			void *param,fitsbtable *fb,int is_pad);

int		fits_bintable_duplicate(fitsbtable *ret,fitsbtable *src,int f);

/* fits-image.c */ /**********************************************************/

/* fits_image_alloc(), fits_image_alloc_2d():
   Allocates an image (see the '*data' elements of the 'fitsimage' structure) 
   and sets the 'sx' and 'sy' elements to the appropriate values.            */
int		fits_image_alloc(fitsimage *fi,int sx,int sy);
int		fits_image_alloc_2d(fitsimage *fi,int sx,int sy);
/* fits_image_alloc_3d():
   Allocates a cubic image (see the 'vdata' and 'allocdata' elements of the 
   'fitsimage' structure) and sets the 'dim' and 'naxis' array elements to the 
   appropriate values. The 'sx' and 'sy' elements are also set and the array 
   'data' points to the first (zeroth) plane in the three dimensional array 
   'vdata'. To access the data, 'vdata' should be casted to (double ***).    */
int		fits_image_alloc_3d(fitsimage *fi,int sx,int sy,int sz);
/* fits_image_alloc_gen():
   Allocates an image with arbitrary dimension 'dim' and with the axis 
   sizes given in the array 'naxis'. The 'sx' and 'sy' elements of the FITS
   structure are also set to 'naxis[0]' and 'naxis[1]' while the array
   'data' points to the first(zeroth) 2 dimensional hyperplane in the
   multidimensional array 'vdata'. If 'dim' is 2, then 'vdata' and 'data'
   has the same value. To access the data, 'vdata' should be casted to
   (double **....**) where the number of asterixes is 'dim'.		     */
int		fits_image_alloc_gen(fitsimage *fi,int dim,int *pnaxis);

double **	fits_image_first_layer(void *data,int dim);
double *	fits_image_first_pixel(void *data,int dim);
int		fits_image_total_pixels(int dim,int *naxis);

void		fits_image_free(fitsimage *fi);

int		fits_image_set_value(fitsimage *fi,double value);
int		fits_image_reset(fitsimage *fi);

/* fits_image_read_line(), fits_image_read_line_cb():
   Reads 'sx' elements from the stream 'f'. The input is treated as a type of 
   'bit' (so, this value should be taken from the header BITPIX before reading)
   and in all cases, it is converted to real (BITPIX=-64, double) values
   and stored in 'line'. The number of bytes read is returned, it is always
   'sx' times |'bit'/8|.                                                     */
int		fits_image_read_line(FILE *f,int sx,int bit,double *line);
int		fits_image_read_line_cb(int (*cb_read)(void *,void *,int),
			void *param,int sx,int bit,double *line);

/* fits_image_get_scale():
   Reads the apropriate scaling and bias value from the headerset 'header'
   using BSCALE and BZERO keywords (if they are present and not ambigous).
   If there are no such headers (and the header is not corrupted) the function
   returns 1.0 and 0.0 as scale and bias value.				     */
int		fits_image_get_scale(fitsheaderset *header,fitsimage *fi,
			double *rbscale,double *rbzero);
/* fits_image_get_params():
   Tries to interpret the headerset 'header as an image (so, NAXIS is 2, and
   has a header of BITPIX, NAXIS1 and NAXIS2), and sets the 'bit', 'sx' and 
   'sy' elments in the FITS structure to the appropriate values. 
   If the header really describes an image, the 'sx', 'sy' and 'bit' parameters
   are set and the call returns 0. Otherwise, these parameters of the FITS
   structure are reset (all of them has a value of 0), and the call returns
   a non-zero value. If in the header of 'img' BSCALE and BZERO keywords are 
   found (see fits_get_scale()), the appropriate values are stroed in in the 
   structure img->read and img->curr also. If there are no such keywords, 
   these parameters are set to 1.0 an 0.0 (see fits_image_get_scale()).	     */
int		fits_image_get_params(fitsheaderset *header,fitsimage *fi);

/* fits_image_rescale():
   Rescales the image in 'fi' ('fi->vdata) using the scaling and bias
   parameters in the structure 'fi'->curr. After the rescaling, 
   'fi'->curr.bscale and 'fi'->curr.bzero are set to 1.0 and 0.0 
   respectively. The call does not affect and use the values in 'fi'->read.  */
int		fits_image_rescale(fitsimage *fi);
/* fits_image_backscale():
   Scales back the image in 'fi' using the scaling and bias parameters
   'scale' and 'zero'. These new scaling parameters are also stored in
   the structure 'fi'->curr while the structure 'fi'->read is not affected.
   The usual call of this function is something like 
   	fits_image_backscale(fi,fi->read.bscale,fi->read.bzero);
   Note the differece between the scalings 'fi'->curr and 'fi'->read. The 
   former is used by these functions, the latter is exported and imported
   by the calls fits_image_get_params() and fits_image_set_params().         */
int		fits_image_backscale(fitsimage *fi,double scale,double zero);

/* fits_image_quantize():
   This function rounds the pixels of the image to 'nquantizebit' significant
   bits. If 'nquantizebit' is zero, pixel values are simply rounded down
   to the nearest integer (see floor()), while if 'nquantizebit' is negative,
   the values are rounded down to the nearest multiple of 2^{-nquantizebit}. */
int		fits_image_quantize(fitsimage *fi,int nquantizebit);

/* fits_image_set_params():
   Exports the 'fi'->bit, 'fi'->dim, 'fi'->naxis, 'fi'->read.bscale and 
   'fi'->read.bzero values into the header elements of the headerset 'header'
   named BITPIX, NAXIS, NAXIS..., BSCALE and BZERO respectively.	     */
int		fits_image_set_params(fitsheaderset *header,fitsimage *fi);

/* fits_image_read(), fits_image_read_cb():
   Tries to read the whole image from the stream 'fr' to the 'fi'->vdata array.
   Before this call, the 'img'->dim, 'img'->naxis[] and 'img'->bit values must 
   be set properly and the 'img'->vdata array must exist (must be allocated).*/
int		fits_image_read(FILE *fr,fitsimage *fi);
int		fits_image_read_cb(int (*cb_read)(void *,void *,int),
			void *param,fitsimage *fi);

/* fits_image_skip(), fits_image_skip_cb():
   Seeks the file pointer of 'fr' after the image. 			     */
int		fits_image_skip(FILE *fr,fitsimage *fi);
int		fits_image_skip_cb(int (*cb_read)(void *,void *,int),
			void *param,fitsimage *fi);

/* fits_image_write_line(), fits_image_write_line_cb():
   Writes the pixel data 'line' with a size of 'sx' to the stream 'fw'. The
   data is converted, according to the value of 'bit'. The function returns
   the number of bytes written to the stream.                                */
int		fits_image_write_line(FILE *fw,int sx,int bit,double *line);
int		fits_image_write_line_cb(int (*cb_write)(void *,void *,int),
			void *param,int sx,int bit,double *line);
/* fits_image_write(), fits_image_write_cb():
   Writes the whole image 'fi' to the stream 'fw'. The function returns
   the number of bytes written to the stream. The function pads the written 
   data with zeros (to be a multiple of 2880 bytes) if 'is_pad' is true.     */
int		fits_image_write(FILE *fw,fitsimage *fi,int is_pad);
int		fits_image_write_cb(int (*cb_write)(void *,void *,int),
			void *param,fitsimage *fi,int is_pad);

/* fits_image_dump():
   A new array of doubles is created with the approporiate size
   and the img->vdata array is copied here, if the FITS image 'fi' contains 
   a valid image (it can have any dimension). If the dumping was successful, 
   the function returns the pointer pointing to the newly allocated array 
   (which can be released with free()), otherwise it returns NULL.	     */
double *	fits_image_dump(fitsimage *fi);

/* fits_image_duplicate():
   Duplicates the FITS image 'src' (and stores the new image in 'ret' which
   can be uninitialized before this call). If 'flag' is true, the whole
   image is copied, otherwise this function fills up the data array (see 
   'vdata' element of the structure) with zeros.			     */
int		fits_image_duplicate(fitsimage *ret,fitsimage *src,int flag);

/* fits_image_bitpix_cname():
   Returns a static string which points to a canonincal C-like variable type
   name which can be assigned to the appropriate bitpix value. Note that
   if bitpix is 8, it returns "byte" (instead of "char") and if bitpix is 32, 
   it returns "long" even in 64-bit architectures.			     */
char *		fits_image_bitpix_cname(int bitpix);

/* fits_image_trim():
   Trims a section from the image. The outer pixels are replaced by the
   value defined by the 'outer' variable. Actually, it works only for 
   2 dimensional images (otherwise, the function does nothing just return 
   a negative value indicating the error).				     */

int		fits_image_trim(fitsimage *i,int x0,int y0,int nx,int ny,double outer);

/* fits-draw.c */ /***********************************************************/

/* fits_image_draw_[pixel,line,circle]():
   Draw pixels, lines and circles to the image data of the FITS image 'img' 
   with the specified value of 'color'. These routines might be useful in 
   special cases or during debugging. These functions return zero if the 
   drawing was successful, otherwise return a non-zero value.		     */
int	fits_image_draw_pixel(fitsimage *img,int x,int y,double col);
int	fits_image_draw_line(fitsimage *img,int x,int y,int dx,int dy,
		double col,int style);
int	fits_image_draw_circle(fitsimage *img,int x,int y,int r,double col);

/* fits_draw_[pixel,line,circle]():
   The same like above, but the target image is the primary image of the
   FITS data 'img' (if any). 						     */
int	fits_draw_pixel(fits *img,int x,int y,double color);
int	fits_draw_line(fits *img,int x1,int y1,int x2,int y2,double c,int sty);
int	fits_draw_circle(fits *img,int x,int y,int r,double color);

/* fits-ui.c */ /*************************************************************/

/* fits_basename():
   Extract the `basename' of the filename 'name'. If 'name' ends like
   [<n>] where <n> is a positive integer, this integer is stored in 
   'retframeno' (if it is not NULL) and a pointer to a truncated name 
   (omitting the trailing [...]) is returned in a static memory area
   (which cannot be released by free() and may be overwritten by further
   calls). If the name does not contain any trailing "[...]" string,
   the function returns the same pointer 'name', and 0 is stored in
   'retframeno' (again, if it is not NULL).
   This function is useful for handling image names like IRAF and ds9 do so. */
char *	fits_basename(char *name,int *retframeno);

/* fits_header_export_command_line():
   Exports the command line argument list described by 'argc','argv'
   (passed to main(), respectively) into the FITS header of 'img' with
   the keyword 'kw'. If the whole command line is longer than ~70 characters,
   the list will be splitted to more card images (all with the same keyword).
   All nasty characters (parsed by the shell, such wildcards, parentheses, ...)
   in the command line argument are escaped or the whole argument is quoted  */
int	fits_header_export_command_line(fits *img,
		char *kw,char *first,char *prefix,
		int argc,char *argv[]);

/* fits_error():
   The fits_error() function returns a string describing the error code
   passed in the argument 'errcode' which is one of the definied error codes
   FITSERR_* (see the begining of this header file).			     */
char *	fits_error(int errcode);

FILE *	fits_file_open(char *name);
FILE *	fits_file_create(char *name);
int	fits_file_close(FILE *f);

/*****************************************************************************/

int	fits_read_extension_table(FILE *fr,fitsextension *xtn);
int	fits_read_extension_image(FILE *fr,fitsextension *xtn);
int	fits_read_extension_bintable(FILE *fr,fitsextension *xtn);

fits *	fits_read_frame_as_extension(FILE *fr,int frameno);
fits *	fits_seek_frame_to_image(FILE *fr,int frameno);
fits *	fits_read_frame_to_image(FILE *fr,int frameno);

int		fits_read_header(FILE *fr,fits *img);
int		fits_get_header_count(fits *img,char *hdr);
int		fits_get_header_id(fits *img,char *hdr,int cnt);
fitsheader *	fits_get_header(fits *img,char *hdr,int cnt);

int		fits_get_header_as_double(fits *img,char *hdr,
		double *ret,int is_ambigous_allowed);

int		fits_get_gain(fits *img,double *gain);

/* fits_set_header_any(), ..., fits_set_header_boolean():
   These functions are wrappers to fits_headerset_set_any(), ... calls 
   and operate on the primary header unit of the FITS data 'img'.            */
fitsheader *	fits_set_header_any(fits *img ,char *,int,char *c);
int		fits_set_header_integer(fits *,char *,int,int    val,char *c);
int		fits_set_header_double (fits *,char *,int,double val,char *c);
int		fits_set_header_string (fits *,char *,int,char * str,char *c);
int		fits_set_header_boolean(fits *,char *,int,int    val,char *c);

int		fits_delete_header(fits *img,char *hdr,int cnt);
int		fits_delete_all_header(fits *img,char *hdr);

void		fits_copy_full_header(fits *im1,fits *im2);

int		fits_alloc_image(fits *img,int sx,int sy);
int		fits_alloc_image_2d(fits *img,int sx,int sy);
int		fits_alloc_image_3d(fits *img,int sx,int sy,int sz);
int		fits_alloc_image_gen(fits *img,int dim,int *naxis);

int		fits_set_image(fits *img,double value);
int		fits_reset_image(fits *img);

int		fits_get_scale(fits *img,double *rbscale,double *rbzero);
int		fits_get_image_params(fits *img);

int		fits_set_image_params(fits *img);

/* fits_set_standard(), fits_set_origin(), fits_set_extend():
   Exports the mandatory FITS header 'SIMPLE' and the optional 'ORIGIN' 
   and 'EXTEND' to the primary header of 'img'. The header 'SIMPLE' always 
   has boolean type and true ('T') value. In fits_set_standard() if comment 
   is NULL, it is set to "FITS standard".				     */
int		fits_set_standard(fits *img,char *comment);
int		fits_set_origin(fits *img,char *origin,char *comment);
int		fits_set_extend(fits *img,int flag,char *comment);

int		fits_rescale(fits *img);
int		fits_backscale(fits *img,double bscale,double bzero);

/* fits_read_image_line(), fits_read_image():
   Wrappers to fits_image_read_line() and fits_image_read(). The latter 
   call here reads the image into the primary image in the FITS data 'img'.  */
int		fits_read_image_line(FILE *fr,int sx,int bit,double *line);
int		fits_read_image(FILE *fr,fits *img);

/* fits_read(), fits_read_cb():
   Reads a complete FITS file from the stream 'fr' (or using the callback 
   'cb_read' with the pameter 'param') and tries to interpret it
   as "intelligent" as possible. If the stream is not a real FITS file,
   the function returns NULL. Otherwise allocates a new FITS structure,
   reads all headers and stores it to the 'hdrs' array and after it tries
   to interpret the header as a header which describes an image. If it
   is really an image, the function set the 'sx', 'sy' and 'bit' values
   and reads the complete image to 'data'. If the FITS data contains 
   extensions, they will also be read. If it is not an image, then
   the data is also read but stored in 'rawdata'.
   If any of these steps fails, the function returns NULL and the stream
   position is undefined.                                                    */
fits * 		fits_read(FILE *fr);
fits *		fits_read_cb(int (*cb_read)(void *,void *,int),void *param);
/* fits_read_raw():
   Reads a complete FITS file from the stream 'fr'. Even if it was a real 
   image, the data is treated as raw data and stored in 'rawdata'. This is a 
   fast method to read a complete FITS file and read the header structure,
   so this call might be useful in pipelines which manipulate just some 
   of the header elements.                                                   */
fits *		fits_read_raw(FILE *fr);

/* fits_write_header():
   Writes the primary header unit of the FITS data 'img' to the stream 'fw'.
   This call is equivalent to fits_headerset_write(fw,&img->header). This 
   function also returns the total number of bytes written to the stream.    */
int		fits_write_header(FILE *fw,fits *img);

/* fits_write_image_line(), fits_write_image():
   Wrappers to fits_image_write_line() and fits_image_write(). The latter 
   call here writes the primary image in the FITS data 'img'.		     */
int		fits_write_image_line(FILE *fw,int sx,int bit,double *line);
int		fits_write_image(FILE *fw,fits *img);

/* fits_write(), fits_write_cb(), fits_mem_write():
   Writes the complete FITS data 'img' (including primary header, primary image
   and all extensions with their headers) to the stream 'fw'. The function 
   returns the number of bytes written to the stream. 			     */
int		fits_write(FILE *fw,fits *img);
int		fits_write_cb(int (*cb_write)(void *,void *,int),
			void *param,fits *img);

int		fits_mem_write(void **buffer,fits *img);

double *	fits_dump_image_raw(fits *img);		

fitsextension	*fits_extension_add(fits *img,int count);
fitsextension	*fits_extension_new(fits *img,int type);

/*****************************************************************************/

/* fits_cb_read(), fits_cb_write():
   Callbacks for reading from and writing to files (FILE *), the parameter
   'param' is the file pointer itself casted to void *.			     */
int		fits_cb_read  (void *param,void *dst,int length);
int		fits_cb_write (void *param,void *src,int length);

/*****************************************************************************/

#endif

/*****************************************************************************/
                                       
