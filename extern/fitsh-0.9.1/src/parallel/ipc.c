/*****************************************************************************/
/* ipc.c								     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Simple parallelization API using basic IPC (inter-process communication)  */
/* calls like fork(), socketpair(), select(), ...			     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* (c) 2006; Pal, A. (apal@szofi.elte.hu)				     */
/*****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <signal.h>
#include <unistd.h>
#include <time.h>
#include <sys/select.h>
#include <sys/time.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <sys/wait.h>

#include "parallel.h"

/*****************************************************************************/

typedef struct child	child;
typedef struct status	status;

struct child
 {	int	pid;
	int	hsock;
	status	*st;
 };

struct status
 { 	child	*chsent;
	int	is_finished,rd;
	void	*pnt;
 };

#define		STOP_SIGNAL		SIGUSR1

/*****************************************************************************/

void child_loop_stop(int param)
{
 exit(0);
}

int child_loop(int size,int (*funct)(const void *param,void *i),const void *param,int hsock)
{
 fd_set	fd; 
 void	*curr;
 int	rd,l,t;

 signal(STOP_SIGNAL,child_loop_stop);

 curr=(void *)malloc(size);
 rd=0;
 while ( 1 )
  {	FD_ZERO(&fd);
	FD_SET(hsock,&fd);
	select(hsock+1,&fd,NULL,NULL,NULL);
	l=size-rd;
	t=read(hsock,(void *)((char *)curr+rd),l);
	if ( t<=0 )	/* stream has been closed */
		break;
	else
	 {	rd+=t;
		if ( rd>=size )
		 {	funct(param,curr);
			(void)write(hsock,curr,size);
			rd=0;
		 }
	 }
  };
 free(curr);
 return(0);
}

int create_child(child *c,int size,int (*funct)(const void *param,void *i),const void *param)
{
 int	pid,socks[2];

 socketpair(AF_LOCAL,SOCK_STREAM,0,socks);

 pid=fork();

 if ( pid<0 )			/* error  */
  	return(-1);

 else if ( pid>0 )		/* parent */
  {	c->pid=pid;
	close(socks[0]);
	c->hsock=socks[1];
	c->st=NULL;
	return(0);
  }
 else
  {	close(socks[1]);	/* child  */
	child_loop(size,funct,param,socks[0]);
	close(socks[0]);
 	close(0);	
	close(1);
	close(2);
	exit(0);
	return(0);
  }	

}

int parallel_ipc(void *base,int nmemb,int size,int (*funct)(const void *param,void *i),const void *param,int flags,int nthread)
{
 int	ret,i,j,maxh,t,l; 
 child	*childs,*c;	/* okay, it's 'children'... */
 int	nchild,nfini;
 status	*sts,*ws;
 fd_set	fd;

 if ( nthread<=1 )
  {	ret=parallel_none(base,nmemb,size,funct,param,flags);
	return(ret);
  }

 nchild=nthread;

 childs=(child  *)malloc(sizeof(child )*nchild);
 sts   =(status *)malloc(sizeof(status)*nmemb);

 for ( i=0 ; i<nchild ; i++ )
  {	c=&childs[i];
	create_child(c,size,funct,param);
  }

 for ( i=0 ; i<nmemb ; i++ )
  {	ws=&sts[i];
	ws->chsent=NULL;
	ws->pnt=(void *)((char *)base+i*size);
	ws->is_finished=0;
	ws->rd=0;
  }

 nfini=0;
 while ( nfini<nmemb ) 
  {	for ( i=0 ; i<nchild ; i++ )
	 {	c=&childs[i];
		if ( c->st==NULL && nfini<nmemb )
		 {	for ( j=0 ; j<nmemb ; j++ )
			 {	ws=&sts[j];
				if ( ws->chsent==NULL )
				 {	(void)write(c->hsock,ws->pnt,size);
					ws->chsent=c;
					ws->rd=0;
					ws->is_finished=0;
					c->st=ws;
					break;
				 }
			 }
		 }
	 }
	FD_ZERO(&fd);
	maxh=0;
	for ( i=0 ; i<nchild ; i++ )
	 {	c=&childs[i];
		FD_SET(c->hsock,&fd);
		if ( c->hsock > maxh )	maxh=c->hsock;
	 }
	select(maxh+1,&fd,NULL,NULL,NULL);
	for ( i=0 ; i<nchild ; i++ )
	 {	c=&childs[i];
		if ( ! FD_ISSET(c->hsock,&fd) || c->st==NULL )	continue;
		l=size-c->st->rd;
		t=read(c->hsock,(void *)((char *)c->st->pnt+c->st->rd),l);
		c->st->rd+=t;
		if ( c->st->rd >= size )
		 {	c->st->is_finished=1;
			c->st=NULL;
			nfini++;
		 }
	 }
  };

 for ( i=0 ; i<nchild ; i++ )
  {	c=&childs[i];
	close(c->hsock);
	kill(c->pid,STOP_SIGNAL);
	waitpid(c->pid,&ret,0);
  }

 free(sts);
 free(childs);

 return(0);
}

/*****************************************************************************/
                                                      
                   
