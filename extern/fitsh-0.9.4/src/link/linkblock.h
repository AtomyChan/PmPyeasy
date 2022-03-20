/*****************************************************************************/
/* linkblock.h								     */
/*****************************************************************************/

#ifndef	__LINKBLOCK_H_INCLUDED
#define	__LINKBLOCK_H_INCLUDED	1

typedef struct linkblock linkblock;

struct linkblock
 {	int		x1,y1;
	int		x2,y2;
	linkblock	*next,*prev;
	int		id,flag;
 };

/* linkblock_connect():
   Connects intersecting blocks of the array 'lblocks'. After this call, 
   the structure elements 'next' and 'prev' are set pointing to the 
   next and previous connected elements (inside this array). Note that the
   elements of the array 'lblocks' are presumably re-ordered; the structure
   members 'id' are set before all sortings to ascending order (begining with 
   zero) but the pointers ('next', 'prev') points inside the sorted array.   */
int	linkblock_connect(linkblock *lblocks,int lnblock);

#endif
                                                                             
