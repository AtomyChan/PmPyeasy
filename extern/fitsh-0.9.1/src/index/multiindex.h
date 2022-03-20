/*****************************************************************************/
/* multiindex.h								     */
/*****************************************************************************/

#ifndef	__MULTIINDEX_H_INCLUDED
#define	__MULTIINDEX_H_INCLUDED	1

/*****************************************************************************/

typedef struct
 {	int	depth;
	int	value;
 } nodeindex;

typedef struct
 {	nodeindex	i;
	int		min;
	int		max;
 } nodeinterval;		/* interval is always [min,max) ! */

typedef struct
 {	int	nbit,size;
	void	*param;
	int	*i_pri,*i_sec;
	int	*ind_x;
	int	**ind_xy;
 } multiindex;

/*****************************************************************************/

int	get_nbit(int x);

int	index_quicksort(int *index,int n,
		int (*compare)(int,int,void *),void *param);

int	check_if_in_interval(void *arr,int p,
		int (*compare)(int,double,void *),double *xl,double *xr);

int	search_index_boundaries(void *arr,int *index,int size,
		int (*compare)(int,double,void *),
		double *rx0,int *rleft ,nodeindex *nleft,
		double *rx1,int *rright,nodeindex *nright);

int	get_nodeintervals(int ileft,int iright,int size,nodeinterval *nis);

int	index_subsort(int *index,int (*compare)(int,int,void *),
		int level,int min,int max,void *arr);

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

int	multiindex_create(void *param,int size,
		int (*compare_x)(int,int,void *),
		int (*compare_y)(int,int,void *),
		multiindex *mi);

int	multiindex_range_query(multiindex *mi,
		int (*compare_value_x)(int,double,void *),
		int (*compare_value_y)(int,double,void *),
		double xl,double xr,double yl,double yr,
		int *ret,int maxret);

int	multiindex_free(multiindex *mi);
int	multiindex_reset(multiindex *mi);

/*****************************************************************************/

#endif

/*****************************************************************************/
                                                                  
