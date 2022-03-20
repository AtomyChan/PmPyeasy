/*****************************************************************************/
/* parallel.h								     */
/*****************************************************************************/

#ifndef	__PARALLEL_H_INCLUDED
#define	__PARALLEL_H_INCLUDED	1

/*****************************************************************************/

#define		PARALLEL_TYPE_NONE	0
#define		PARALLEL_TYPE_IPC	1

/*****************************************************************************/

#define		PARALLEL_DEFAULT	0x00
#define		PARALLEL_BREAK		0x01

typedef struct
 {	int		nthread;
 } par_ipc;

typedef struct
 { 	int		_dummy;
 } par_none;

typedef union
 {	par_ipc		ipc;
	par_none	none;
 } par_general;

typedef struct
 { 	int		type;
	par_general	p;
 } parallel;

/*****************************************************************************/

/* parallel_ipc():
   Parallelizes the evaluation of the function 'funct' on the array 'base'.
   The array 'base' has 'nmemb' element, each with a size of 'size'
   (similarly like qsort()). The parallelization is local, using 
   threads created by fork(). The callback function 'funct' has 
   two arguments: the optional parameter 'param' and the array element
   pointer 'dpnt'. The argument 'flags' is for fine-tuning the algorithm,
   while the argument 'nthread' specifies the number of threads to use
   (in the actual implementation, 'nthread' new processes are created, 
   the parent is related to the control of the whole thing...)		     */
int parallel_ipc(void *base,int nmemb,int size,
	int (*funct)(const void *param,void *dpnt),const void *param,
	int flags,int nthread);

/* parallel_none():
   This function has the same calling syntax like above (except for the 
   argument 'nthread'), does the same but without parallelization. This
   function might be useful for testing the behaviour of 'funct()' in
   a single-process environment.					     */
int parallel_none(void *base,int nmemb,int size,
	int (*funct)(const void *param,void *dpnt),const void *param,
	int flags);

/*
int parallel_cb_ipc(int nmemb,int size,
*/

/*****************************************************************************/

int parallel_general(void *base,int nmemb,int size,
	int (*funct)(const void *param,void *dpnt),const void *param,
	int flags,parallel *pg);

/*****************************************************************************/

#endif

/*****************************************************************************/
                                                       
                     
