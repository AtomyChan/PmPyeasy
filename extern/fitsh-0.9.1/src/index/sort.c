/*****************************************************************************/
/* sort.c 								     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Standalone library for sorting an indexed array.			     */
/* (c) 2004, Pal, A. (apal@szofi.elte.hu).				     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* See function prototypes and the usage of the functions in sort.h.	     */
/*****************************************************************************/

#include <stdio.h>
#include <stdlib.h>

#include "sort.h"

/*****************************************************************************/

static int index_qsort_local(int *index,int l,int r,
	int (*compare)(int,int,void *),void *param)
{
 int	i,j,tmp,ch;

 ch=index[(l+r)/2];
 i=l-1,j=r+1;
 while ( 1 )
  {	do { i++; } while ( index[i] != ch && compare(index[i],ch,param)<0 ) ;
	do { j--; } while ( index[j] != ch && compare(index[j],ch,param)>0 ) ;
	if ( i<j )
		tmp=index[j],index[j]=index[i],index[i]=tmp;
	else if ( i==j )
	 {	i++;
		break;
	 }
	else	break;
  };
 if ( l<j )	index_qsort_local(index,l,j,compare,param);
 if ( i<r )	index_qsort_local(index,i,r,compare,param);
 return(0);
}
int index_qsort_old(int *index,int n,int (*compare)(int,int,void *),void *param)
{
 index_qsort_local(index,0,n-1,compare,param);
 return(0);
}

/*****************************************************************************/
 
typedef struct
 {	int	lo;
	int	hi;
 } stacknode;

#define		MAX_THRESH	4

#define		PUSH(low,high)	((void)((top->lo=(low)),(top->hi=(high)),top++))
#define		POP(low,high)	((void)(top--,(low=top->lo),(high=top->hi)))
                                           
int index_qsort(int *index,int n,int (*compare)(int,int,void *),void *param)
{
 int	t;

 if ( n<=0 )	return(0);

 if ( n>MAX_THRESH )
  {	stacknode	stack[64],*top;
	int		lo=0,hi=n-1,mid,left,right;
	top=&stack[1];

	while (	top>stack )
	 {	mid=lo+((hi-lo)>>1);
		if ( compare(index[mid],index[lo],param)<0 )
			t=index[mid],index[mid]=index[lo],index[lo]=t;
		if ( compare(index[hi],index[mid],param)<0 )
			t=index[hi],index[hi]=index[mid],index[mid]=t;
		else	goto	skip_compare;
		if ( compare(index[mid],index[lo],param)<0 )
			t=index[mid],index[mid]=index[lo],index[lo]=t;
		skip_compare:
		left=lo+1;
		right=hi-1;
		do
		 {	while ( compare(index[left],index[mid],param)<0 )
				left++;
			while ( compare(index[mid],index[right],param)<0 )
				right--;
			if ( left<right )
			 {	t=index[left],index[left]=index[right],index[right]=t;
				if ( mid==left )
					mid=right;
				else if ( mid==right )
					mid=left;
				left++;
				right--;
			 }
			else if ( left==right )
			 {	left++;
				right--;
				break;
			 }
		 } while ( left<=right );

		if ( right-lo <= MAX_THRESH )
		 {	if ( hi-left <= MAX_THRESH )
				POP(lo,hi);
			else
				lo=left;
		 }
		else if ( hi-left <= MAX_THRESH )
			hi=right;
		else if ( right-lo > hi-left )
		 {	PUSH(lo,right);
			lo=left;
		 }
		else
		 {	PUSH(left,hi);
			hi=right;
		 }
	 }
  }

 do	
  {	int	end=n-1;
	int	tmp=0,run,trav;
	int	thresh=(end>MAX_THRESH?MAX_THRESH:end);

	for ( run=tmp+1 ; run<=thresh ; run++ )
	 {	if ( compare(index[run],index[tmp],param)<0 )
			tmp=run;
	 }
	if ( tmp != 0 )
		t=index[0],index[0]=index[tmp],index[tmp]=t;
	run=1;
	while ( (run+=1) <= end )
	 {	tmp=run-1;
		while ( compare(index[run],index[tmp],param)<0 )
			tmp--;
		tmp++;
		if ( tmp != run )
		 {	trav=run+1;
			while ( (--trav)>=run )
			 {	int	c=index[trav];
				int	hi,lo;
				for ( hi=lo=trav ; (lo-=1)>=tmp ; hi=lo )
					index[hi]=index[lo];
				index[hi]=c;
			 }
		 }
	 }
  } while(0);

 return(0);
}

/*****************************************************************************/
                                               
