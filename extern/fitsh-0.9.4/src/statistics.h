/*****************************************************************************/
/* statistics.h 								     */
/*****************************************************************************/

#ifndef	__STATISTICS_H_INCLUDED
#define	__STATISTICS_H_INCLUDED 1

/*****************************************************************************/

/* median():
   Sorts the array 'data' of 'n' real values into ascending order and returns
   the median of the dataset.						     */
double	median(double *data,int n);

/* mean():
   Mean of the array 'data'.						     */
double	mean(double *data,int n);

/* truncated_mean():
   Sorts the array 'data' of 'n' real values into ascending order and returns
   the truncated mean of the dataset by rejecting the lower and upper 't' points */
double	truncated_mean(double *data,int n,int t);


/*****************************************************************************/

#endif

/*****************************************************************************/
                                
