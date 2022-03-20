/*****************************************************************************/
/* tensor.h 								     */
/*****************************************************************************/

#ifndef	__TENSOR_H_INCLUDED
#define	__TENSOR_H_INCLUDED	1

/* tensor_alloc_arr():
   Allocates memory for a tensor with rank 'rank', the number of the elements 
   in each dimension have to be specified in the elements of 'dims' ('rank' 
   integers). If the allocation is successful, the function returns the pointer
   pointing to the newly allocated tensor, otherwise (or if 'rank'<=0) it 
   returns NULL. If a tensor of a given 'type' wants to be allocated, set 
   'typesize' to sizeof('type') and after the call, recast the pointer returned
   to ('type' **...**), where the indirection-level (the number of the 
   asterices) is always 'rank'.						     */
void *	tensor_allor_arr(int typesize,int rank,int *dims);

/* tensor_alloc():
   Does the same what tensor_alloc_arr() does, except that the number of 
   elements in each dimension should be specified as integers after the 
   argument 'rank' (also, 'rank' integers).				     */
void *	tensor_alloc(int typesize,int rank,...);

/* tensor_free():
   Releases the tensor 'tensor' allocated with tensor_alloc[_arg](). In the
   current implementation, it calls directly free(), however, it might be 
   changed in the further versions. Currently, this function does always 
   return zero.								     */
int	tensor_free(void *tensor);

#define	tensor_alloc_1d(type,a)		tensor_alloc(sizeof(type),1,a)
#define	tensor_alloc_2d(type,a,b)	tensor_alloc(sizeof(type),2,a,b)
#define	tensor_alloc_3d(type,a,b,c)	tensor_alloc(sizeof(type),3,a,b,c)
#define	tensor_alloc_4d(type,a,b,c,d)	tensor_alloc(sizeof(type),4,a,b,c,d)

#endif
                         
