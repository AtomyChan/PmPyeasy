/*****************************************************************************/
/* cache.c								     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Standalone library for chunked I/O caching on large files.		     */ 
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* (c) 2006; Pal, A. (apal@szofi.elte.hu)				     */
/*****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <string.h>
#include <errno.h>

#ifndef	HOST_WIN32
#define	USE_MMAP
#endif

#ifdef	USE_MMAP
#include <sys/mman.h>
#endif

#include "cache.h"

/*****************************************************************************/

static int get_gcd(int a,int b)	/* greatest common divisor */
{
 int	c;
 if ( a<=0 || b<=0 )	return(-1);
 while ( (c=a%b)>0 )	a=b,b=c;
 return(b);
}

static int get_lcm(int a,int b)	/* least common multiple */
{
 if ( a<=0 || b<=0 )	return(-1);
 return(a*b/get_gcd(a,b));
}

/*****************************************************************************/

int cache_getpagesize(void)
{
#ifdef	HOST_WIN32
 return(4096);
#else
 return(getpagesize());
#endif
}

int cache_blocksize(int recordsize)
{
 int	pagesize,blocksize;

 pagesize=cache_getpagesize();
 blocksize=get_lcm(recordsize,pagesize);

 return(blocksize);
}

int cache_tcblock(int recordsize,off_t recordcount,int multip)
{
 int	pagesize,chpagesize,chnrecord,tcblock;

 if ( multip<=0 )	multip=1;

 pagesize=cache_getpagesize();
 chpagesize=get_lcm(recordsize,pagesize)*multip;
 chnrecord=chpagesize/recordsize;
 tcblock=(int)((recordcount+(off_t)chnrecord-(off_t)1)/(off_t)chnrecord);
 return(tcblock);
}

int cache_init(cache *ch,int recordsize,off_t recordcount,int multip,
	int nchblock,int fh,int is_write)
{
 int		i;
 cacheblock	*cb;

 if ( multip<=0 )	multip=1;
 if ( nchblock<=0 )	nchblock=1;

 ch->fh=fh;
 ch->is_write=is_write;

 ch->recordsize=recordsize;
 ch->recordcount=recordcount;

 ch->pagesize=cache_getpagesize();
 ch->chpagesize=get_lcm(ch->recordsize,ch->pagesize)*multip;
 ch->chnrecord=ch->chpagesize/ch->recordsize;

 ch->tcblock=(int)((recordcount+(off_t)ch->chnrecord-(off_t)1)/(off_t)ch->chnrecord);
 ch->nn_lookup=(short *)malloc(sizeof(short)*ch->tcblock);
 for ( i=0 ; i<ch->tcblock ; i++ )
  {	ch->nn_lookup[i]=-1;		}
 
 ch->nchblock=nchblock;
 ch->chblocks=(cacheblock *)malloc(sizeof(cacheblock)*ch->nchblock);
 ch->nusedblock=0;
 ch->nmapped=0;
 for ( i=0 ; i<ch->nchblock ; i++ ) 
  {	cb=&ch->chblocks[i];
	cb->cachearray=NULL;
	cb->size=0;
	cb->nnch=-1;
	cb->prev=NULL;
	cb->next=NULL;
  }

 ch->recent=NULL;
 ch->oldest=NULL;

 return(0);
}

/*****************************************************************************/

static cacheblock *cache_get_block(cache *ch,int nnch,cacheblock **rchrej)
{
 cacheblock	*cb,*rc,*rj;

 if ( ch->recent==NULL )
  {	cb=&ch->chblocks[0];
	ch->recent=cb;
	ch->oldest=cb;
	cb->prev=NULL;
	cb->next=NULL;
	ch->nusedblock=1;
	cb->nnch=nnch;
	ch->nn_lookup[nnch]=0;
	if ( rchrej != NULL )	*rchrej=NULL;
	return(cb);
  }
 
 if ( ch->recent->nnch==nnch )
  {	if ( rchrej != NULL )	*rchrej=NULL;
	return(ch->recent);
  }
 if ( ch->nn_lookup[nnch]>=0 )
  {	rc=&ch->chblocks[ch->nn_lookup[nnch]];		}
 else
	rc=NULL;
	
 if ( rc != NULL )
  {	if ( rc->prev != NULL )	rc->prev->next=rc->next;
	else			ch->oldest=rc->next;
	if ( rc->next==NULL )
	 {	fprintf(stderr,"cache_get_block():127: rc->next == NULL\n");
		fprintf(stderr,"nnch=%d\n",nnch);
		fprintf(stderr,"ch->recordcount=%ld\n",(long)ch->recordcount);
		fprintf(stderr,"ch->recent->nnch=%d\n",ch->recent->nnch);
		fprintf(stderr,"ch->nusedblock=%d,ch->nchblock=%d\n",ch->nusedblock,ch->nchblock);
		exit(1);
	 }
	rc->next->prev=rc->prev;
	ch->recent->next=rc;
	rc->prev=ch->recent;
	rc->next=NULL;
	ch->recent=rc;
	rc->nnch=nnch;
	if ( rchrej != NULL )	*rchrej=NULL;
	return(rc);
  } 
 else if ( ch->nusedblock<ch->nchblock )
  {	rc=&ch->chblocks[ch->nusedblock];
	ch->nn_lookup[nnch]=ch->nusedblock;
	ch->nusedblock++;
	ch->recent->next=rc;
	rc->prev=ch->recent;
	rc->next=NULL;
	ch->recent=rc;
	rc->nnch=nnch;
	if ( rchrej != NULL )	*rchrej=NULL;
	return(rc);
  }
 else
  {	rj=ch->oldest;
	ch->nn_lookup[rj->nnch]=-1;
	rj->next->prev=NULL;
	rc=ch->oldest;
	ch->oldest=rj->next;
	ch->recent->next=rc;
	rc->prev=ch->recent;
	rc->next=NULL;
	ch->recent=rc;
	rc->nnch=nnch;
	ch->nn_lookup[nnch]=(ch->recent)-(ch->chblocks);
	if ( rchrej != NULL )	*rchrej=rj;
	return(rc);
  }
 
}

static cacheblock * cache_do_mapping(cache *ch,off_t recnum)
{
 int		nnch;
 cacheblock	*cb,*cj;
 off_t		offset,wf;
 int		size;

 nnch=recnum/(off_t)ch->chnrecord;
 cb=cache_get_block(ch,nnch,&cj);

 /* this block should be released: */
 if ( cj != NULL && cj->cachearray != NULL )	
  {
#ifdef USE_MMAP
	munmap((void *)cj->cachearray,(size_t)cj->size);
#else
	if ( ch->is_write )
	 {	lseek(ch->fh,cj->offset,SEEK_SET);
		write(ch->fh,cj->cachearray,cj->size);
	 }
	free(cj->cachearray);
#endif
	ch->nmapped--;
	cj->cachearray=NULL;
	cj->size=0;
  }

 if ( cb->cachearray==NULL )
  {	offset=(off_t)ch->chpagesize*(off_t)nnch;
	wf=(off_t)ch->recordsize*(off_t)ch->recordcount-offset;
	if ( wf>(off_t)ch->chpagesize )	size=ch->chpagesize;
	else				size=(int)wf;
	cb->size=size;
	cb->offset=offset;
	if ( ch->is_write )	
	 {	
#ifdef	USE_MMAP
		cb->cachearray=(char *)mmap(NULL,(size_t)cb->size,
			PROT_READ|PROT_WRITE,MAP_SHARED,ch->fh,offset);
		if ( (void *)cb->cachearray==(void *)(-1) ) /* mmap() failed */
		 {	fprintf(stderr,"[rw] offset=%d %p %d nnch=%d chnrecord=%d fh=%d error='%s'\n",
			(int)offset,(void *)cb->cachearray,(int)cb->size,nnch,ch->chnrecord,ch->fh,strerror(errno));
			return(NULL);
		 }
#else
		cb->cachearray=(char *)malloc((size_t)cb->size);
		lseek(ch->fh,(off_t)offset,SEEK_SET);
		read(ch->fh,cb->cachearray,(size_t)cb->size);
#endif
		ch->nmapped++;
	 }
	else
	 {	
#ifdef	USE_MMAP
		cb->cachearray=(char *)mmap(NULL,(size_t)cb->size,
			PROT_READ,MAP_SHARED,ch->fh,offset);
		if ( (void *)cb->cachearray==(void *)(-1) ) /* mmap() failed */
		 {	fprintf(stderr,"[ro] offset=%d %p %d nnch=%d chnrecord=%d fh=%d error='%s'\n",
			(int)offset,(void *)cb->cachearray,(int)cb->size,nnch,ch->chnrecord,ch->fh,strerror(errno));
			return(NULL);
		 }
#else
		cb->cachearray=(char *)malloc((size_t)cb->size);
		lseek(ch->fh,(off_t)offset,SEEK_SET);
		read(ch->fh,cb->cachearray,(size_t)cb->size);
#endif
		ch->nmapped++;
	 }
  }

 return(cb);
}

void  *cache_read_record(cache *ch,off_t recnum)
{
 int		sbrc;
 cacheblock	*cb;
 void		*record;

 cb=cache_do_mapping(ch,recnum);
 if ( cb==NULL )	return(NULL);	/* cache_do_mapping() failed */

 sbrc=(int)(recnum%(off_t)ch->chnrecord);
 record=(void *)(cb->cachearray+sbrc*ch->recordsize);

 return(record);
}

int cache_write_record(cache *ch,off_t recnum,void *record)
{
 int		sbrc;
 cacheblock	*cb;

 cb=cache_do_mapping(ch,recnum);
 if ( cb==NULL )	return(-1);	/* cache_do_mapping() failed */

 sbrc=(int)(recnum%(off_t)ch->chnrecord);
 memcpy((void *)(cb->cachearray+sbrc*ch->recordsize),record,ch->recordsize);

 return(0);
}

/*****************************************************************************/

static int cache_do_unmapping(cache *ch,cacheblock *cb)
{
 if ( ch==NULL || cb==NULL )
	return(-1);
 if ( cb->cachearray==NULL || cb->nnch<0 )
	return(1);

#ifdef	USE_MMAP
 munmap((void *)cb->cachearray,(size_t)cb->size);
#else
 if ( ch->is_write )
  {	lseek(ch->fh,cb->offset,SEEK_SET);
	write(ch->fh,cb->cachearray,cb->size);
  }
 free(cb->cachearray);
#endif

 return(0);
}

static int cache_free(cache *ch)
{
 if ( ch->chblocks != NULL )
	free(ch->chblocks);
 if ( ch->nn_lookup != NULL )
	free(ch->nn_lookup);
 return(0);
}
int cache_finalize(cache *ch)
{
 int		i;
 cacheblock	*cb;
 for ( i=0 ; i<ch->nchblock ; i++ )
  {	cb=&ch->chblocks[i];
	cache_do_unmapping(ch,cb);
  }
 cache_free(ch);
 return(0);
}

/******************************************************************************
int fprint_chain_forward(FILE *fw,cache *ch)
{
 cacheblock	*cb;
 if ( ch->oldest==NULL )	fprintf(fw,"<NULL>\n");
 else  
  {	for ( cb=ch->oldest ; cb != NULL ; cb=cb->next )
	 {	fprintf(fw,"<%d> ",cb->nnch);			}
	fprintf(fw,"\n");
  }
 return(0);
}
int fprint_chain_backward(FILE *fw,cache *ch)
{
 cacheblock	*cb;
 if ( ch->recent==NULL )	fprintf(fw,"<NULL>\n");
 else  
  {	for ( cb=ch->recent ; cb != NULL ; cb=cb->prev )
	 {	fprintf(fw,"<%d> ",cb->nnch);			}
	fprintf(fw,"\n");
  }
 return(0);
}
int main(int argc,char *argv[])
{
 cache		ch;
 char		buff[256];
 int		nnch;
 cacheblock	*cb;
 cache_init(&ch,0,256,65536,5);
 while ( 1 ) 
  {	fprint_chain_forward(stdout,&ch);
	fprint_chain_backward(stdout,&ch);
	nnch=-1;
	while ( nnch<0 )
	 {	fgets(buff,255,stdin);
		sscanf(buff,"%d",&nnch);
	 };
	cb=cache_get_block(&ch,nnch,NULL);
	cb->nnch=nnch;
  };
}
******************************************************************************/

typedef struct
 {	off_t	lo;
	off_t	hi;
 } stacknode;

#define	MAX_THRESH		4
#define	PUSH(low,high)		((void)((top->lo=(low)),(top->hi=(high)),top++))
#define	POP(low,high)		((void)(top--,(low=top->lo),(high=top->hi)))
#define	COMPARE(_p1,_p2,_v1,_v2)	(_v1=cache_read_record(ch,(_p1)),\
					 _v2=cache_read_record(ch,(_p2)),\
					 compare(_v1,_v2))
#define	SWAP(_p1,_p2,_v1,_v2)		(memcpy(t1,_v1,recordsize),\
					 memcpy(t2,_v2,recordsize),\
					 cache_write_record(ch,(_p1),t2),\
					 cache_write_record(ch,(_p2),t1))

int cache_sort_block(cache *ch,off_t o,off_t n,
	int (*compare)(const void *,const void *))
{
 void		*vl,*vr,*vm,*t1,*t2;
 stacknode	stack[64],*top;
 off_t		lo,hi,mid,left,right,end,tmp,run,trav,thresh;
 int		recordsize;

 if ( n<=0 )	return(0);

 recordsize=ch->recordsize;

 t1=(void *)malloc(recordsize);
 t2=(void *)malloc(recordsize);

 lo=o;
 hi=o+n-1;

 if ( n>MAX_THRESH )
  {	top=&stack[1];
	while ( top>stack )
	 {	mid=lo+((hi-lo)>>1);
		if ( COMPARE(mid,lo,vm,vl)<0 )
			SWAP(mid,lo,vm,vl);
		if ( COMPARE(hi,mid,vr,vm)<0 )
			SWAP(hi,mid,vr,vm);
		else	goto	skip_compare;
		if ( COMPARE(mid,lo,vm,vl)<0 )
			SWAP(mid,lo,vm,vl);
		skip_compare:

		left=lo+1;
		right=hi-1;
		do
		 {	while ( COMPARE(left,mid,vl,vm)<0 )
				left++;
			while ( COMPARE(mid,right,vm,vr)<0 )
				right--;
			if ( left<right )
			 {	COMPARE(left,right,vl,vr);
				SWAP(left,right,vl,vr);
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
 end=n-1;
 tmp=0;
 thresh=(end<MAX_THRESH?end:MAX_THRESH);
 for ( run=tmp+1 ; run<=thresh ; run++ )
  {	if ( COMPARE(o+run,o+tmp,vl,vr)<0 )
		tmp=run;
  }
 if ( tmp != 0 )
  {	COMPARE(o+0,o+tmp,vl,vr);
	SWAP(o+0,o+tmp,vl,vr);
  }
 run=o+1;
 while ( (run+=1) <= end )
  {	tmp=run-1;
	while ( COMPARE(o+run,o+tmp,vl,vr)<0 )
		tmp--;
	tmp++;
	if ( tmp != run )
	 {	trav=run+1;
		while ( (--trav)>=run )
		 {	off_t	hi,lo;
			memcpy(t1,cache_read_record(ch,o+trav),recordsize);
			for ( hi=lo=trav ; (lo-=1)>=tmp ; hi=lo )
			 {	memcpy(t2,cache_read_record(ch,o+lo),recordsize);
				cache_write_record(ch,o+hi,t2);
			 }
			cache_write_record(ch,o+hi,t1);
		 }
	 }
  }

 free(t2);
 free(t1);

 return(0);
}

int cache_sort(cache *ch,int (*compare)(const void *,const void *))
{
 return(cache_sort_block(ch,0,ch->recordcount,compare));
}

/*****************************************************************************/
                                             
