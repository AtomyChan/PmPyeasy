/*****************************************************************************/
/* random.c								     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* (C) 2011; Pal, A. (apal@szofi.net); based on the code available in 	     */
/* the book: Numerical Recipes: the art of scientific computing (3rd ed.)    */
/*****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

#include <random/random.h>

static struct	random_state	state;
static int	is_initialized=0;

#define		R2P32	2.3283064365386962890625E-10

/*****************************************************************************/

int random_state_seed(struct random_state *state,int j)
{
 state->v  = 2244614371U;
 state->w1 =  521288629U;
 state->w2 =  362436069U;
 state->u  =  (uint32_t)j ^ state->v;
 random_state_uint32(state);
 state->v  = state->u;
 random_state_uint32(state);
 return(0);
}

uint32_t random_state_uint32(struct random_state *state)
{
 uint32_t	x,y;

 state->u = state->u * 2891336453U + 1640531513U;
 state->v ^= state->v >> 13;
 state->v ^= state->v << 17;
 state->v ^= state->v >>  5;
 state->w1 = 33378 * (state->w1 & 0xffff) + (state->w1 >> 16);
 state->w2 = 57225 * (state->w2 & 0xffff) + (state->w2 >> 16);
 x = state->u ^ (state->u << 9); 
 x ^= x >> 17; 
 x ^= x <<  6;
 y = state->w1 ^ (state->w1 << 17); 
 y ^= y >> 15; 
 y ^= y << 5;
 return((x+state->v)^(y+state->w2));
}

double random_state_double(struct random_state *state)
{
 uint32_t	r1,r2;
 r1=random_state_uint32(state);
 r2=random_state_uint32(state);
 return(R2P32*((double)r1+R2P32*(double)r2));
}

/*****************************************************************************/

int random_seed(int j)
{
 random_state_seed(&state,j);
 is_initialized=1;
 return(0);
}

uint32_t random_uint32(void)
{ 
 if ( ! is_initialized )
  {	random_seed(0);
	is_initialized=1;
  }
 return(random_state_uint32(&state));
}

double random_double(void)
{
 if ( ! is_initialized )
  {	random_seed(0);
	is_initialized=1;
  }
 return(random_state_double(&state));
}

/*****************************************************************************/
     
