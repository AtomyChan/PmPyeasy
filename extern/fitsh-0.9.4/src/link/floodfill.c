/*****************************************************************************/
/* floodfill.c 								     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Standalone library for filling areas.				     */
/* (c) 2004, Pal, A. (apal@szofi.elte.hu).	 			     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* See function prototypes and the usage of the functions in floodfill.h     */
/*****************************************************************************/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#include "floodfill.h"

/*****************************************************************************/

typedef struct
 {	int	x1,x2;
	int	y,dy;
 } linesegment;

#define MAXDEPTH 64

#define PUSH(XL,XR,Y,DY) \
	if ( sp<stack+MAXDEPTH ) \
	 {	sp->x1=XL;sp->x2=XR;sp->y=Y;sp->dy=DY; sp++;	}

#define POP(XL,XR,Y,DY) \
	 {	sp--;XL=sp->x1;XR=sp->x2;Y=sp->y+(DY=sp->dy);	}

int floodfill_old(int x,int y,int (*getp)(void *,int,int),void (*putp)(void *,int,int),void *param)
{
 int		left,x1,x2,dy,skip;
 linesegment	stack[MAXDEPTH],*sp;

 sp=stack;

 if ( getp(param,x,y) )	return(0);

 PUSH(x,x,y,1);
 PUSH(x,x,y+1,-1);

 left=0;
 while ( sp>stack )
  {	POP(x1,x2,y,dy);

	for ( x=x1 ; ! getp(param,x,y) ; x-- )
		putp(param,x,y);

	if( x >= x1 )	skip=1;
	else		skip=0;
	
	if ( ! skip )
	 {	left=x+1;
		if ( left<x1 )
		 {	PUSH(y,left,x1-1,-dy);		}
		x=x1+1;
	 }

	do
	 {	if ( ! skip )
		 {	for( ; ! getp(param,x,y) ; x++ )
				putp(param,x,y);
			PUSH(left,x-1,y,dy);
			if ( x>x2+1 )
			 {	PUSH(x2+1,x-1,y,-dy);	}
		 }
		skip=0;
		for( x++ ; x <= x2 && getp(param,x,y) ; x++ ) ;
		left=x;
	 } while ( x<=x2 );
  }
 return(0);
}

typedef struct 
 {	int	x;
	int	y;
 } point;

#define		STACKBLOCK		64

int floodfill(int x,int y,int (*getp)(void *,int,int),void (*putp)(void *,int,int),void *param)
{
 point	stack_static[STACKBLOCK],*stack_dynamic,*stack;
 int	bl,bh,px,py,is_dec,sp,mxsp;

 if ( getp(param,x,y) )	return(0);

 stack=stack_static;
 stack_dynamic=NULL;
 sp=0;mxsp=STACKBLOCK;

 bl=bh=0;
 is_dec=1;

 while ( 1 )
  {	if ( is_dec )
	 {	while ( ! getp(param,x,y) )	x--;	}
	putp(param,x,y);
	if ( ! getp(param,x,y-1) )
	 {	if ( bl==0 )
		 {	bl=1;
			stack[sp].x=x,stack[sp].y=y-1,sp++;
			if ( sp==mxsp )
			 {	mxsp+=STACKBLOCK;
				if ( sp==STACKBLOCK )
				 {	stack_dynamic=(point *)malloc(sizeof(point)*mxsp);
					memcpy(stack_dynamic,stack_static,sizeof(point)*sp);
					stack=stack_dynamic;
				 }
				else
				 {	stack_dynamic=(point *)realloc(stack_dynamic,sizeof(point)*mxsp);
					stack=stack_dynamic;
				 }
			 }
		 }
	 }
	else	bl=0;
	if ( ! getp(param,x,y+1) )
	 {	if ( bh==0 )
		 {	bh=1;
			stack[sp].x=x,stack[sp].y=y+1,sp++;
			if ( sp==mxsp )
			 {	mxsp+=STACKBLOCK;
				if ( sp==STACKBLOCK )
				 {	stack_dynamic=(point *)malloc(sizeof(point)*mxsp);
					memcpy(stack_dynamic,stack_static,sizeof(point)*sp);
					stack=stack_dynamic;
				 }
				else
				 {	stack_dynamic=(point *)realloc(stack_dynamic,sizeof(point)*mxsp);
					stack=stack_dynamic;
				 }
			 }
		 }
	 }
	else	bh=0;
	if ( ! getp(param,x+1,y) )
	 {	x++;
		is_dec=0;
	 }
	else
	 {	if ( sp==0 )	break;
		sp--,px=stack[sp].x,py=stack[sp].y;
		if ( y<py || ( y==py && x<px ) )	bh=0,x=px,y=py;
		else					bl=0,x=px,y=py;
		is_dec=1;
	 }
  }

 if ( stack_dynamic != NULL )	free(stack_dynamic);

 return(0);
}

/*****************************************************************************/
          
