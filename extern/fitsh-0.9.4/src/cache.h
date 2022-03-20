/*****************************************************************************/
/* cache.c								     */
/*****************************************************************************/

#ifndef	__CACHE_H_INCLUDED
#define	__CACHE_H_INCLUDED	1

#include <sys/types.h>		/* definition of type off_t */

/*****************************************************************************/

typedef struct	cacheblock	cacheblock;

struct cacheblock
 {	char		*cachearray;
	off_t		offset;
	int		size;
	int		nnch;
	cacheblock	*prev;
	cacheblock	*next;
 };

typedef struct
 {	int		fh;		/* file handle			     */
	int		is_write;	/* writeable cache		     */

	int		recordsize;
	off_t		recordcount;	

	int		pagesize;	/* see getpagesize()		     */
	int		chpagesize;	/* LCM of pagesize and recordsize    */
	int		chnrecord;	/* chpagesize / recordsize	     */
					/* = num of record in a cache block  */
	int		tcblock;	/* total cache blocks		     */
	
	cacheblock 	*chblocks;

	int		nchblock,	/* should be less than 32k	     */
			nusedblock,
			nmapped;

	short		*nn_lookup;	/* lookup table for cache blocks,    */
					/* this array has 'tcblock' elements */

	cacheblock	*recent;
	cacheblock	*oldest;
 } cache;

/*****************************************************************************/

#define		CACHE_READONLY		(0)	/* ro - cache */
#define		CACHE_READ_AND_WRITE	(1)	/* rw - cache */

/* cache_blocksize():
   Returns the size of the smallest cache block, which is the least common 
   multiple of recordsize  and the architecture's preferred pagesize 
   (reported by the getpagesize() system call).				     */
int	cache_blocksize(int recordsize);

/* cache_tcblock():
   calculates the required value of 'tcblock' which is enough to hold 
   all data in the cache at the same time (if it is passed to cache_init()
   as the value of 'nchblock' (see cache_init() for more details).	     */
int	cache_tcblock(int recordsize,off_t recordcount,int multip);

/* cache_init():
   Initializes the cache 'ch' on the opened file descriptor 'fh'. 
   The number of 'recordcount' objects (each with the size of 'recordsize')
   are assumed to be stored in the file. If we want to write into the 
   file (using caching), set 'is_write' to a non-zero value (or
   CACHE_READ_AND_WRITE), otherwise set it to zero (or CACHE_READONLY).
   The caching will be done in 'nchblock' cache blocks, each block will have
   the size of 'multip' times the least common multiple of recordsize
   and the architecture's preferred pagesize (reported by the getpagesize()
   system call). The total count of cache blocks which will be used can be
   figured out using cache_tcblock(). If 'nchblock' is less than or equal to
   this value, all data can be cached into memory. Note that 'nchblock' cannot
   exceed 32k-1 (32767). If all data wants to be cached, increase 'multip'
   to have a tcblock count less than 32k.
     Note also that the maximum number of 'recordcount' is also architecture
   and/or compile-time dependent. If "off_t" is a 32-bit integer (sizeof(off_t) 
   is 4), then the maximum number of objects is 2G, if "off_t" is a 64-bit 
   integer (sizeof(off_t) is 8), then the maximum number of objects is 8E. The 
   implementation of "off_t" depends on the architecture and can be tuned
   during compilation time using -D_FILE_OFFSET_BITS={32|64}. Be careful,
   because of even on modern unixes and 64-bit architectures(!), the default
   "off_t" type can either be a 32-bit integer. This library not uses the 
   type of "long long", therefore it is fully ANSI-compatible and pedantic.  */
int	cache_init(cache *ch,int recordsize,off_t recordcount,
		int multip,int nchblock,int fh,int is_write);

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/* cache_read_record():
   Returns a pointer to the 'recnum'th record using the cache 'ch'. Simple.
   The only important thing is that the subsequent calls of cache_read_record()
   may overwrite the internal buffers and the pointers previously returned 
   by this call will have an unexpected contain. Because of the implementation,
   these contains are assured to be preserved if the total number of subsequent
   calls (including cache_write_record()'s also) is less than or equal to 
   'nchblock' (see cache_init() above).					     */
void  *	cache_read_record(cache *ch,off_t recnum);

/* cache_write_record():
   Stores the record 'record' in the cache 'ch' at the 'recnum'th record.    */
int	cache_write_record(cache *ch,off_t recnum,void *record);

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/* cache_sort_block(): 
   Sorts the first 'n' record in the cache 'ch', from the offset 'o'. 
   Note that o+n cannot be larger than ch->recordcount.                      */
int	cache_sort_block(cache *ch,off_t o,off_t n,
	int (*compare)(const void *,const void *));

/* cache_sort(): 
   Sorts the contents of the cache, using the function 'compare' for 
   comparison.								     */
int	cache_sort(cache *ch,int (*compare)(const void *,const void *));


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/* cache_finalize():
   Flushes and releases all cache data in 'ch'. The underlying file won't be
   closed by this call.							     */
int	cache_finalize(cache *ch);

/*****************************************************************************/

#endif

/*****************************************************************************/
               
