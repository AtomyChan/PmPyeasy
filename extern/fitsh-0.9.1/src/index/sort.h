/*****************************************************************************/
/* sort.h 								     */
/*****************************************************************************/

#ifndef	__SORT_H_INCLUDED
#define	__SORT_H_INCLUDED	1

/*****************************************************************************/

/* index_qsort():
   Sorts the array 'index' to satisfy the relation 
	   compare(index[i],index[i+1],param) <= 0     (for all i=0,... , n-2).
   The function index_qsort() uses the quick sort algorihm enhanced with
   insertion sort on shorter blocks of the array. The optional argument 
   'param' is always passed to the comparison function 'compare' at
   all invocations.							     */
int index_qsort(int *index,int n,
	int (*compare)(int index1,int index2,void *param),void *param);

/* index_qsort_old():
   Previous implementation. Obsoleted by index_qsort().			     */
int index_qsort_old(int *index,int n,
	int (*compare)(int index1,int index2,void *param),void *param);

/*****************************************************************************/

#endif

/*****************************************************************************/
                        
