/*****************************************************************************/
/* random.h								     */
/*****************************************************************************/

#ifndef	__RANDOM_H_INCLUDED
#define	__RANDOM_H_INCLUDED	1

#include <stdint.h>

/*****************************************************************************/

struct random_state
 {	uint32_t	u,v;
	uint32_t	w1,w2;
 };

int		random_state_seed(struct random_state *state,int j);
uint32_t 	random_state_uint32(struct random_state *state);
double 		random_state_double(struct random_state *state);

int		random_seed(int j);
uint32_t 	random_uint32(void);
double 		random_double(void);

/*****************************************************************************/

#endif

/*****************************************************************************/
                                                           
