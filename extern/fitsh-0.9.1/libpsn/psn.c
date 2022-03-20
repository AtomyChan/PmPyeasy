/*****************************************************************************/
/* psn.c 								     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Standalone library for symbolic mathematical calculations based on	     */
/* PSN (Polish Suffix Notation) sequences.				     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Version: 1.0rc3, Last modified: 2009.05.24.				     */
/* (c) 1996, 2002-2003, 2004-2005, 2007; Pal, A. (apal@szofi.elte.hu). 	     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* See function prototypes and the usage of the functions in psn.h	     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* This program is free software; you can redistribute it and/or modify	     */
/* it under the terms of the GNU General Public License as published by      */
/* the Free Software Foundation; either version 2 of the License, or	     */
/* (at your option) any later version.					     */
/*									     */
/* This library is distributed in the hope that it will be useful,	     */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of	     */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the	     */
/* GNU General Public License for more details.				     */
/*									     */
/* You should have received a copy of the GNU General Public License	     */
/* along with this program; if not, write to the Free Software		     */
/* Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111, USA.    */
/*****************************************************************************/

#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#include <psn/psn.h>

#define		psncpy(a,b)	memcpy((void *)(a),(void *)(b),sizeof(psnterm))

#define		FRAGSIZE_TERM		32
#define		FRAGSIZE_CONS		8
#define		PSNSTACKBLOCK		64
#define		BLOCKSIZE		128

/*****************************************************************************/

int	psnerrno=0;

/*****************************************************************************/

typedef struct
 {	int	offset;
	int	length;
 } istack;

typedef struct
 {	istack	orig;
	istack	diff;
 } dstack;

/*****************************************************************************/

static psn * psn_alloc(void)
{
 psn *ret;
 ret=(psn *)malloc((unsigned)sizeof(psn));
 if ( ret==NULL )	return(NULL);
 ret->terms=(psnterm *)malloc(sizeof(psnterm)*FRAGSIZE_TERM);
 if ( ret->terms==NULL )
  {	free(ret);return(NULL);		}
 ret->ntermalloc=FRAGSIZE_TERM;
 ret->nterm=0;
 ret->terms[0].type=T_END;

 ret->cons=NULL;
 ret->ncon=0;
 ret->nconalloc=0;
 ret->nseq=0;

 return(ret);
}

void psn_free(psn *seq)
{
 if ( seq==NULL )		return;
 if ( seq->terms != NULL )	free(seq->terms);
 if ( seq->cons != NULL )	free(seq->cons);
 free(seq);
}

static void psn_memory_compact(psn *seq)
{
 int k;
 k=seq->nterm+1;
 k=((k+FRAGSIZE_TERM-1)/FRAGSIZE_TERM)*FRAGSIZE_TERM;
 if ( k < seq->ntermalloc )
  {	seq->terms=(psnterm *)realloc(seq->terms,k*sizeof(psnterm));
	seq->ntermalloc=k;
  }
}
void psn_append(psn *seq,psnterm *terms,int n)
{
 int k;
 k=seq->nterm+1+n;
 if ( k>seq->ntermalloc )
  {	k=seq->ntermalloc=((k+FRAGSIZE_TERM-1)/FRAGSIZE_TERM)*FRAGSIZE_TERM;
	seq->terms=(psnterm *)realloc(seq->terms,k*sizeof(psnterm));
  }
 memcpy(seq->terms+seq->nterm,terms,n*sizeof(psnterm));
 seq->nterm+=n;
 seq->terms[seq->nterm].type=T_END;
 seq->terms[seq->nterm].major=0;
}
static void psn_append_inseq(psn *seq,int termpos,int n)
{
 int k;
 k=seq->nterm+1+n;
 if ( k>seq->ntermalloc )
  {	k=seq->ntermalloc=((k+FRAGSIZE_TERM-1)/FRAGSIZE_TERM)*FRAGSIZE_TERM;
	seq->terms=(psnterm *)realloc(seq->terms,k*sizeof(psnterm));
  }
 memmove(seq->terms+seq->nterm,seq->terms+termpos,n*sizeof(psnterm));
 seq->nterm+=n;
 seq->terms[seq->nterm].type=T_END;
 seq->terms[seq->nterm].major=0;
}
void psn_append_term(psn *seq,int type,int major,int minor,int cache)
{
 psnterm term;

 term.type =type,
 term.major=major,
 term.minor=minor,
 term.cache=cache;
 psn_append(seq,&term,1);
}
static void psn_remove(psn *seq,int st,int num)
{
 int l;

 l=seq->nterm;
 if ( st+num>l || num==0 )	return;
 memmove(seq->terms+st,seq->terms+st+num,(l-st-num+1)*sizeof(psnterm));
 seq->nterm-=num;
 seq->terms[seq->nterm].type=T_END;
 seq->terms[seq->nterm].major=0;
}
static void psn_insert(psn *seq,int st,int num)
{
 int k;
 k=seq->nterm+1+num;
 if ( k>seq->ntermalloc )
  {	k=seq->ntermalloc=((k+FRAGSIZE_TERM-1)/FRAGSIZE_TERM)*FRAGSIZE_TERM;
	seq->terms=(psnterm *)realloc(seq->terms,k*sizeof(psnterm));
  }
 memmove(seq->terms+st+num,seq->terms+st,(seq->nterm-st+1)*sizeof(psnterm));
 seq->nterm+=num;
 seq->terms[seq->nterm].type =T_END;
 seq->terms[seq->nterm].major=0;
}

static void psn_swap_sequences(psn *p1,psn *p2)
{
 psnterm *w;
 int	 l;

 w=p1->terms     ,p1->terms     =p2->terms     ,p2->terms     =w;
 l=p1->nterm     ,p1->nterm     =p2->nterm     ,p2->nterm     =l;
 l=p1->ntermalloc,p1->ntermalloc=p2->ntermalloc,p2->ntermalloc=l;
}

static void psn_cons_alloc(psn *seq)
{
 seq->cons=(double *)malloc(sizeof(double)*FRAGSIZE_CONS);
 seq->ncon=0,seq->nconalloc=FRAGSIZE_CONS;
}
int psn_cons_append(psn *seq,double con)
{
 int	ret;
 if ( seq->cons==NULL )	psn_cons_alloc(seq);

 if ( seq->ncon+1>seq->nconalloc )
  {	seq->nconalloc+=FRAGSIZE_CONS;
	seq->cons=(double *)realloc(seq->cons,sizeof(double)*seq->nconalloc);
  }
 ret=seq->ncon;
 seq->cons[ret]=con;
 seq->ncon++;
 return(ret);
}
static void psn_cons_append_more(psn *seq,double *cons,int n)
{
 while ( n )
  {	psn_cons_append(seq,*cons);
	cons++,n--;
  };
}

/*****************************************************************************/

static psnprop* psn_prop_search(psnprop *pl,int major)
{
 for ( ; pl->major ; pl++ )
  {	if ( pl->major == major )	return(pl);	}
 return(NULL);
}

static int psn_get_precedency(psnterm *ps1,psnterm *ps2,psnprop *plist)
{
 psnprop *pl1,*pl2;
 int	p1,p2;

 for ( pl1=pl2=NULL ; plist->major != 0 && (pl1==NULL||pl2==NULL) ; plist++ )
  {	if ( ps1->major==plist->major )	pl1=plist;
	if ( ps2->major==plist->major )	pl2=plist;
  }

 p1=p2=0;
 if ( pl1 != NULL )	p1=pl1->precedency;
 if ( pl2 != NULL )	p2=pl2->precedency;
 if ( ps1->type==T_FN )	p1=PSN_MAX_PREC;
 if ( ps2->type==T_FN )	p2=PSN_MAX_PREC;

 if ( p1==0 || p2==0 )	return(0);

 if ( p1==p2 )
  {	if ( p1==PSN_MAX_PREC )	return(1);
	else if ( pl1==pl2 )	return(pl1->associativity);
	else			return(0);
  }
 else	return(p1-p2);
}

psn * psn_conv(psn *chseq,psnprop *plist)
{
 int		sp,plevel;
 psn		*ret;
 psnterm	*in,*stack,stack_loc[PSNSTACKBLOCK],*stack_dyn;
 int		stacklen;

 in=chseq->terms;

 ret=psn_alloc();
 ret->nseq=1;

 if ( chseq->cons != NULL )
	psn_cons_append_more(ret,chseq->cons,chseq->ncon);

 stacklen=PSNSTACKBLOCK,sp=0;
 stack_dyn=NULL;
 stack=stack_loc;

 plevel=0;

 while ( in->type )
  {	if ( in->type==T_CONST || in->type==T_SCONST || in->type==T_VAR )
	 {	psn_append(ret,in,1);		 }

	else if ( in->type==T_PAROPEN )
	 {	psncpy(stack+sp,in);sp++,plevel++;	}

	else if ( ( in->type==T_PARCLOSE || in->type==T_SEP ) && plevel>0 )
	 {	sp--;
		if ( in->type==T_PARCLOSE ) plevel--;
		while ( stack[sp].type != T_PAROPEN )
		 {	if ( stack[sp].type==T_FN )	stack[sp].type=T_OP;
			psn_append(ret,stack+sp,1);
			sp--;
		 };
		if ( in->type==T_SEP )
		 {	psncpy(stack+sp,in);
			stack[sp].type=T_PAROPEN;sp++;
		 }
	 }
	else if ( in->type==T_SEP && plevel==0 )
	 {	while ( sp )
		 {	sp--;
			if ( stack[sp].type==T_FN )	stack[sp].type=T_OP;
			psn_append(ret,stack+sp,1);
		 };
		ret->nseq++;
	 }

	else if ( in->type==T_OP || in->type==T_FN )
	 {	
		if ( sp==0 )
		 {	psncpy(stack+sp,in);sp++;	}

		else if ( (stack+sp-1)->type==T_PAROPEN )
		 {	psncpy(stack+sp,in);sp++;	}

		else if ( psn_get_precedency(in,stack+sp-1,plist) > 0 )  
		 {	psncpy(stack+sp,in);sp++;	}

		else
		 {	while ( psn_get_precedency(stack+sp-1,in,plist)>=0 )
			 {	sp--;

				if ( stack[sp].type==T_FN )
					stack[sp].type=T_OP;

				psn_append(ret,stack+sp,1);

				if ( sp==0 || stack[sp-1].type==T_PAROPEN )
					break;
			 };
			psncpy(stack+sp,in);sp++;
		 }
	 }
	else	
	 {	psnerrno=PEINVALID;
		if ( stack_dyn != NULL )	free(stack_dyn);
		return(NULL);
	 }

	if ( sp>=stacklen )
	 {	if ( stacklen==PSNSTACKBLOCK )
		 {	stacklen+=PSNSTACKBLOCK;
			stack_dyn=(psnterm *)malloc(sizeof(psnterm)*stacklen);
			memcpy(stack_dyn,stack,sizeof(psnterm)*sp);
			stack=stack_dyn;
		 }
		else
		 {	stacklen+=PSNSTACKBLOCK;
			stack=(psnterm *)realloc(stack,sizeof(psnterm)*stacklen);
			stack_dyn=stack;
		 }	
	 }

	in++;
  };

 while ( sp > 0 )
  {	sp--;

	if ( stack[sp].type==T_FN )
		stack[sp].type=T_OP;

	psn_append(ret,stack+sp,1);
  };

 if ( stack_dyn != NULL )	free(stack_dyn);

 psnerrno=PEOK;
 return(ret);
}

/*****************************************************************************/

static int is_cchar_start(int c)
{
 if ( ('a'<=c && c<='z')||('A'<=c && c<='Z')|| c=='_' )
	return(1);
 else	
	return(0);
}

static int is_cchar(int c)
{
 if ( ('a'<=c && c<='z')||('A'<=c && c<='Z')|| c=='_' ||( '0'<=c && c<='9') )
	return(1);
 else	
	return(0);
}

static char opchars[]="!#$%&*+-/:<=>@^|~";

static int is_opchar(int c)
{
 char	*ch;
 for ( ch=opchars ; *ch ; ch++ )
  {	if ( c==(*ch) )	return(1);	}
 return(0);
}
static int is_numchar_start(int c)
{
 if ( ( '0'<=c && c<='9' ) || c=='.' )	return(1);
 else					return(0);
}
static int is_numchar(int c)
{
 if ( ( '0'<=c && c<='9' ) || c=='.' || c=='e' || c=='E' )	return(1);
 else								return(0);
}

static int psn_get_symbol_from_symtable	(char *wp,int wl,psnsym **symtable,
	   psnsym *found,int max)
{
 int	i,nf;

 for ( i=0,nf=0 ; *symtable != NULL && nf<max ; )
  {	if ( (*symtable)[i].type==0 || (*symtable)[i].name==NULL )
	 {	i=0;
		symtable++;
		continue;
	 }
	if ( strlen((*symtable)[i].name)==wl && 
	memcmp((*symtable)[i].name,wp,wl)==0 )
	 {	found[nf].type =(*symtable)[i].type ,
		found[nf].major=(*symtable)[i].major,
		found[nf].minor=(*symtable)[i].minor;
		found[nf].name =NULL;
		nf++;
	 }
	i++;
  }
 return(nf);
}

psn * psn_conv_string(char *str,psnsym **symtable)
{
 int	i,nc,wl,nf;
 char	*wp;
 double	ind;
 short	ins;
 psnsym	found[3];

 psn	*ret;

 ret=psn_alloc();
 ret->nseq=1;

 nc=0;	

 for ( ; *str ; str++ )
  {	if ( *str==32 || *str==10 || *str==13 || *str==9 )	continue;

	else if ( *str=='(' )	psn_append_term(ret,T_PAROPEN ,0,0,0);
	else if ( *str==')' )	psn_append_term(ret,T_PARCLOSE,0,0,0);
	else if ( *str=='{' )	psn_append_term(ret,T_PAROPEN ,1,0,0);
	else if ( *str=='}' )	psn_append_term(ret,T_PARCLOSE,1,0,0);
	else if ( *str=='[' )	psn_append_term(ret,T_PAROPEN ,2,0,0);
	else if ( *str==']' )	psn_append_term(ret,T_PARCLOSE,2,0,0);
	else if ( *str==',' )	psn_append_term(ret,T_SEP,     0,0,0);
	else if ( *str==';' )	psn_append_term(ret,T_SEP,     1,0,0);

	else if ( is_cchar_start(*str) )
	 {	wl=1,wp=str,str++;
		while ( is_cchar(*str) )	str++,wl++;
		str--;
		nf=psn_get_symbol_from_symtable(wp,wl,symtable,found,1);
		if ( nf==0 )
		 {	psnerrno=PENOTFOUND;
			return(NULL);
		 }
		psn_append_term(ret,found[0].type,found[0].major,
			found[0].minor,0);
	 }
	else if ( is_opchar(*str) )
	 {	int	fprf,finf,fsff,k;
		wl=1,wp=str,str++;
		while ( is_opchar(*str) )	str++,wl++;
		str--;
		nf=psn_get_symbol_from_symtable(wp,wl,symtable,found,3);
		if ( nf==0 )
		 {	psnerrno=PENOTFOUND;
			return(NULL);
		 }
		fprf=finf=fsff=-1;
		for ( k=0 ; k<nf ; k++ )
		 {	     if ( found[k].minor<0  )	fprf=found[k].major;
			else if ( found[k].minor==0 )	finf=found[k].major;
			else				fsff=found[k].major;
		 }
		psn_append_term(ret,T_OP,fprf,finf,fsff);
	 }
	else if ( is_numchar_start(*str) )
	 {	wl=1,wp=str,str++;
		while ( is_numchar(*str) || ( ( *str=='+' || *str=='-' ) &&
		( *(str-1)=='e' || *(str-1)=='E' ) ) )
			str++,wl++;
		str--;

		sscanf(wp,"%lg",&ind);ins=(short)ind;
		if ( (double)ins == ind )
			psn_append_term(ret,T_SCONST,ins,0,0);
		else
		 {	for ( i=0 ; i<ret->ncon ; i++ )
			 {	if ( ret->cons[i]==ind )	break;	}
			if ( i<ret->ncon && ret->ncon>0 )
			 	psn_append_term(ret,T_CONST,i,0,0);
			else
			 {	psn_append_term(ret,T_CONST,nc,0,0);
				psn_cons_append(ret,ind);
				nc++;
			 }
		 }
	 }
	else
	 {	psnerrno=PEINVALID;
		return(NULL);		/* illegal character */
	 }
  };

 psnerrno=PEOK;
 return(ret);
}

/*****************************************************************************/

static int psn_init_parentheses(psn *pseq)
{
 int	 k,wk,type,nseq;
 psnterm *seq;
 int	 *stack,stack_local[3*PSNSTACKBLOCK],*stack_dyn,stacklen,sp; 
 
 seq=pseq->terms;
 stacklen=3*PSNSTACKBLOCK;
 stack=stack_local,stack_dyn=NULL;
 sp=0;

 nseq=1;

 for ( k=0 ; seq[k].type ; k++ )
  {	type=seq[k].type;
	if ( type==T_PAROPEN )
	 {	stack[sp+0]=1,
		stack[sp+1]=k;
		stack[sp+2]=seq[k].major;
		sp+=3;
	 }
	else if ( type==T_SEP )
	 {	if ( sp>0 )
		 {	if ( stack[sp-1] != seq[k].major )
				return(PEPARSE);
			stack[sp-3]++;
		 }
		else	
			nseq++;
	 }
	else if ( type==T_PARCLOSE )
	 {	if ( sp==0 )
			return(PEPARSE);
	 	sp-=3;
		if ( stack[sp+2] != seq[k].major )
			return(PEPARSE);
		wk=stack[sp+1];
		if ( k>wk+1 )	seq[wk].minor=seq[k].minor=stack[sp+0];
		else		seq[wk].minor=seq[k].minor=0;
	 }
	if ( sp>=stacklen-3 )
	 {	if ( stacklen==3*PSNSTACKBLOCK )
		 {	stacklen+=3*PSNSTACKBLOCK;
			stack_dyn=(int *)malloc(sizeof(int)*stacklen);
			memcpy(stack_dyn,stack,sizeof(int)*sp);
			stack=stack_dyn;
		 }
		else
		 {	stacklen+=3*PSNSTACKBLOCK;
			stack=(int *)realloc(stack_dyn,sizeof(int)*stacklen);
			stack_dyn=stack;
		 }
	 }
  }

 pseq->nseq=nseq;

 if ( stack_dyn != NULL )	free(stack_dyn);

 if ( sp>0 )	return(PEPARSE);
 else		return(0);
}

static int psn_init_function_arguments(psn *pseq)
{
 int	 k,narg;
 psnterm *seq;

 seq=pseq->terms;

 for ( k=0 ; seq[k].type ; k++ )
  {	if ( seq[k].type==T_FN && seq[k+1].type==T_PAROPEN )
	 {	narg=seq[k+1].minor;
		if ( seq[k].minor < 0 )
			seq[k].minor=narg;
		else if ( seq[k].minor != narg )
			return(PEPARSE);
		else
			seq[k].minor=narg;
	 }
	else if ( seq[k].type==T_FN )
		seq[k].minor=0;
  }
 return(0);
}

static int psn_init_operators(psn *pseq)
{
 int	k,tprev,tfoll,iprf,iinf,isff;
 int	canbe_prefix,canbe_suffix;
 psnterm *seq;

 seq=pseq->terms;

 for ( k=0,tprev=0 ; seq[k].type ; k++ )
  {	tfoll=seq[k+1].type;
	if ( seq[k].type==T_OP )
	 {	iprf=seq[k].major,
		iinf=seq[k].minor,
		isff=seq[k].cache;

		if ( tprev==0 || tprev==T_PAROPEN  || tprev==T_OP ||
		tprev==T_FN || tprev==T_SEP )
			canbe_prefix=1;
		else 
			canbe_prefix=0;

		if ( tfoll==0 || tfoll==T_PARCLOSE || tfoll==T_OP || 
		tfoll==T_SEP )
			canbe_suffix=1;
		else
			canbe_suffix=0;

		if ( canbe_prefix && (!canbe_suffix) )
			seq[k].major=iprf,seq[k].minor=1;

		else if ( canbe_suffix && (!canbe_prefix) )
			seq[k].major=isff,seq[k].minor=1;

		else if ( (!canbe_prefix) && (!canbe_suffix) )
			seq[k].major=iinf,seq[k].minor=2;

		else
			return(PEPARSE);

		if ( seq[k].major < 0 )
			return(PEPARSE);

		else if ( seq[k].major==0 )	/* Remove identity... */
			psn_remove(pseq,k,1);

	 }
	tprev=seq[k].type;
  };
 return(0);
}

int psn_cache_clear(psn *seq)
{
 psnterm *term;

 for ( term=seq->terms ; term->type != 0 ; term++ )
  {	if ( term->type==T_OP || term->type==T_FN )
		term->cache=-1;
  }

 return(0);
}

int psn_init(psn *pseq,psnprop *plist)
{
 int	k;

 k=psn_init_parentheses(pseq); 
 if ( k )	return(psnerrno=k);
 k=psn_init_function_arguments(pseq);
 if ( k )	return(psnerrno=k);
 k=psn_init_operators(pseq);
 if ( k )	return(psnerrno=k);

 psn_cache_clear(pseq);
 psnerrno=PEOK;
 return(0);
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

int psn_cache_search(int major,psnfunct *flist)
{
 int	  k;

 for ( k=0 ; flist[k].major != 0 ; k++ )
  {	if ( major == flist[k].major )
		return(k);
  }
 return(-1);
}

int psn_cache_init(psn *pseq,psnfunct *flist)
{
 psnterm  *seq;
 int	  k;

 for ( seq=pseq->terms ; seq->type != 0 ; seq++ )
  {	if ( seq->type==T_OP || seq->type==T_FN )
	 {	k=psn_cache_search(seq->major,flist);
		if ( k<0 )
			return(psnerrno=PENOTFOUND);
		seq->cache=k;
	 }
  }

 return(psnerrno=0);
}

/*****************************************************************************/

static char * psn_symtable_search(psnsym **symtable,int major,int t1,int t2)
{
 int i;
 if ( symtable==NULL )
	return(NULL);
 while ( *symtable != NULL )
  {	for ( i=0 ; (*symtable)[i].type != 0 ; i++ )
	 {	if ( (  (*symtable)[i].type==t1 || (*symtable)[i].type==t2 ) &&
		(*symtable)[i].major==major )
			return((*symtable)[i].name);
	 }
	symtable++;
  };
 return(NULL);
}
void psn_fprint(FILE *f,psn *seq,psnsym **symtable)
{
 psnterm *psn;
 char *name;

 for ( psn=seq->terms ; psn->type ; psn++ )
  {	
	if ( psn->type==T_PAROPEN )
		fprintf(f,"(/%d/",psn->minor);
	else if ( psn->type==T_PARCLOSE )
		fprintf(f,")");

	else if ( psn->type==T_SEP )
		fprintf(f,",");

	else if ( psn->type==T_CONST )
		fprintf(f,"%g",seq->cons[psn->major]);

	else if ( psn->type==T_SCONST )
		fprintf(f,"%g",(double)(psn->major));

	else if ( psn->type==T_OP || psn->type==T_FN )	
	 {	name=psn_symtable_search(symtable,psn->major,T_OP,T_FN  );
		if ( name != NULL )	fprintf(f,"[%s]/%d/",name,psn->minor);
		else			fprintf(f,"???");
	 }
	else if ( psn->type==T_VAR )	
	 {	name=psn_symtable_search(symtable,psn->major,T_VAR,T_VAR);
		if ( name != NULL )	fprintf(f,"%s",name);
		else			fprintf(f,"???");
	 }
	else if ( psn->type==T_STACKVAR )	fprintf(f,"(%d)",psn->major);
	fprintf(f," ");
  }
}

/*****************************************************************************/

int psn_get_nstack(psn *pseq,int *ropt)
{
 int	 sp,nopt,narg;
 psnterm *seq;

 seq=pseq->terms;
 nopt=0;

 for ( sp=0 ; seq->type ; seq++ )
  {	switch( seq->type )
	 {   case T_CONST: case T_SCONST: case T_VAR: 
		sp++;
		break;
	     case T_STACKVAR:
		if ( seq->major+1>nopt )	nopt=seq->major+1;
		sp++;
		break;
	     case T_OP: case T_FN:
		narg=seq->minor;
		if ( sp<narg )
	 	 {	psnerrno=PEPARSE;
			return(-1);
		 }
		sp+=1-narg;
		break;
	 }
  }
 psnerrno=0;
 if ( ropt != NULL )	*ropt=nopt;
 return(sp);
}

int psn_get_length(psn *seq)
{
 return(seq->nterm);
}
int psn_get_nseq(psn *seq)
{
 return(seq->nseq);
}

int psn_test(psn *pseq)
{
 int	nopt,n;

 n=psn_get_nstack(pseq,&nopt);
 if ( n<0 )	return(n);

 if ( n==pseq->nseq+nopt )	return(psnerrno=0);
 else				return(psnerrno=PEPARSE);
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

int psn_double_calc(psn *pseq,psnfunct *flist,double *result,double *vars)
{
 double		stack_automatic[PSNSTACKBLOCK],*stack_dynamic,*stack;
 int		i,nopt,narg,sp;
 psnterm	*seq;
 int		stacklen;

 stack=stack_automatic;
 stack_dynamic=NULL;
 stacklen=PSNSTACKBLOCK;

 nopt=0;sp=0;
 for ( seq=pseq->terms ; seq->type ; seq++ )
  {	switch( seq->type )
	 {   case T_CONST:
		stack[sp]=pseq->cons[seq->major],sp++;
	  	break;
	     case T_VAR:
		stack[sp]=vars[seq->major],sp++;
		break;
	     case T_SCONST:
		stack[sp]=(double)(seq->major),sp++;
		break;
	     case T_STACKVAR:
		stack[sp]=stack[seq->major],sp++;
		if ( seq->major+1>nopt )	nopt=seq->major+1;
		break;
	     case T_OP: case T_FN:
		narg=seq->minor;
		i=seq->cache;
		if ( i<0 )
		 {	i=psn_cache_search(seq->major,flist);
			if ( i<0 )
			 {	if ( stack_dynamic != NULL )
					free(stack_dynamic);
				return(psnerrno=PENOTFOUND);
			 }
			seq->cache=i;
		 }
		if ( flist[i].funct(stack+sp) )
		 {	if ( stack_dynamic != NULL )
				free(stack_dynamic);
			return(psnerrno=PENUMERICAL);
		 }
		sp+=1-narg;
		break;
	 }
	if ( sp>=stacklen )
	 {  if ( stacklen==PSNSTACKBLOCK )
	     {	int	i;
		stacklen+=PSNSTACKBLOCK;
		stack_dynamic=(double *)malloc(sizeof(double)*stacklen);
		for ( i=0 ; i<PSNSTACKBLOCK ; i++ )
		 {	stack_dynamic[i]=stack_automatic[i];	}
	     }
	    else
	     {	stacklen+=PSNSTACKBLOCK;
		stack_dynamic=(double *)realloc(stack_dynamic,
		sizeof(double)*stacklen);
	     }
	    stack=stack_dynamic;
	 }
  }

 for ( i=0 ; i<sp-nopt ; i++ )
  {	result[i]=stack[nopt+i];		}

 if ( stack_dynamic != NULL )
	free(stack_dynamic);

 return(psnerrno=0);
}

/*****************************************************************************/

static int psn_diff_check_var(psnterm *terms,int n,int var)
{
 for ( ; n>0 ; terms++,n-- )
  {	if ( terms->type==T_VAR && terms->major==var )
		return(1);
  }
 return(0);
}

static psn * psn_diff_int(psn *pseq,psnprop *plist,psndiff *rules,
		int var,int flags)
{
 psn	 *ncc,*dcc;
 psnterm *seq;
 psnprop *pr;
 int	 i,j,k,pn,pd,sp,*pc,stacklen,in,major;
 short	 *dr;

 for ( i=0 ; i<pseq->nterm ; i++ )
  {	if ( pseq->terms[i].type==T_STACKVAR )
	 {	psnerrno=PEOPTIMIZED;
		return(NULL);
	 }
  }

 dcc=psn_alloc();
 if ( dcc==NULL )
  {	psnerrno=PEALLOC;
	return(NULL);
  }
 ncc=psn_alloc();
 if ( ncc==NULL )
  {	psnerrno=PEALLOC;
	psn_free(dcc);
	return(NULL);
  }
 ncc->nseq=pseq->nseq;

 if ( pseq->cons != NULL )
	psn_cons_append_more(dcc,pseq->cons,pseq->ncon);

 stacklen=PSNSTACKBLOCK;
 pc=(int *)malloc(sizeof(int)*4*stacklen);
 sp=0;pc[0]=pc[1]=pc[2]=pc[3]=0;
 
 for ( seq=pseq->terms ; seq->type ; seq++ )
  {	
	switch( seq->type )
	 {   case T_CONST: case T_SCONST: case T_VAR:

		if ( seq->type==T_VAR && var==seq->major )  major=1;
		else					    major=0;
		psn_append_term(dcc,T_SCONST,major,0,0);
		pd=pc[4*sp];pd++;

		psn_append(ncc,seq,1);
		pn=pc[4*sp+2];pn++;

		pc[4*sp+1]=pc[4*sp+3]=1;

		sp++;

		pc[4*sp]=pd,pc[4*sp+2]=pn;
		pc[4*sp+1]=pc[4*sp+3]=0;
		break;

	     case T_OP: case T_FN:

		pr=psn_prop_search(plist,seq->major);
		if ( pr==NULL )
		 {	psn_free(ncc);psn_free(dcc);
			psnerrno=PENOTFOUND;
			return(NULL);
		 }
		in=pr->argnum;		/* number of operands	*/

		for ( i=0 ; rules[i].major != 0 ; i++ )
		 {	if ( rules[i].major==seq->major )  break;  }

		if ( rules[i].major==0 )
	 	 {	if ( ! flags )
			 {	psn_free(ncc);psn_free(dcc);
				psnerrno=PENONDIFF;
				return(NULL);
			 }
			else
				dr=NULL;
		 }
		else
		 	dr=rules[i].oplist;	/* get differental rule */

		pd=pc[4*sp];

		if ( dr != NULL )
		 {	i=0;
			for ( ; *dr ; dr++ )
			 {	if ( *dr>=0 )
			 	 {	pr=psn_prop_search(plist,*dr);
					if ( pr==NULL )
					 {	psn_free(ncc);psn_free(dcc);
						psnerrno=PENOTFOUND;
						return(NULL);
					 }
					psn_append_term(dcc,T_OP,*dr,pr->argnum,0);
					i++;
				 }
				else if ( -SS_REF+1<=*dr && *dr<=-1 )
				 {	j=(-(*dr))- 0;k=pc[4*(sp-j)+1];
					psn_append_inseq(dcc,pc[4*(sp-j)+0],k);
					i+=k;
				 }
				else if ( -2*SS_REF+1<=*dr && *dr<=-SS_REF-1 )
				 {	j=(-(*dr))-SS_REF;k=pc[4*(sp-j)+3];
					psn_append(dcc,ncc->terms+pc[4*(sp-j)+2],k);
					i+=k;
				 }
				else
				 {	j=(-(*dr))-2*SS_REF;
					psn_append_term(dcc,T_SCONST,j,0,0);
					i++;
				 }
			 }
		 }
		else
		 {	for ( i=0 ; i<in ; i++ )
			 {	k=pc[4*(sp-in+i)+3];
				if ( psn_diff_check_var(ncc->terms+pc[4*(sp-in+i)+2],k,var) )
				 {	psn_free(ncc);psn_free(dcc);
					psnerrno=PENONDIFF;
					return(NULL);
				 }
			 }
			psn_append_term(dcc,T_SCONST,0,0,0);
			i=1;
		 }

		pn=pc[4*sp+2];
		psn_append(ncc,seq,1);

		sp-=in;	

		psn_remove(dcc,pc[4*sp],pd-pc[4*sp]);
		pc[4*sp+1]=i;
		
		pc[4*sp+3]++;
		for ( i=1 ; i<in ; i++ )
		 {	pc[4*sp+3]+=pc[4*(sp+i)+3];	}

		sp++;
		pc[4*sp+0]=pc[4*(sp-1)+0]+pc[4*(sp-1)+1];
		pc[4*sp+2]=pc[4*(sp-1)+2]+pc[4*(sp-1)+3];
		pc[4*sp+1]=pc[4*sp+3]=0;
		
		break;
	 }
	if ( sp>=stacklen-2 )
	 {	stacklen+=PSNSTACKBLOCK;
		pc=(int *)realloc(pc,sizeof(int)*4*stacklen);
	 }
  }

 psn_free(ncc);

 psn_memory_compact(dcc);

 dcc->nseq=pseq->nseq;

 psn_cache_clear(dcc);

 free(pc);

 psnerrno=PEOK;
 return(dcc);
}

psn * psn_diff(psn *pseq,psnprop *plist,psndiff *rules,int var)
{
 psn	*diff;
 diff=psn_diff_int(pseq,plist,rules,var,0);
 return(diff);
}

/*****************************************************************************/

static int psn_int_cmp(psn *seq,int p1,int p2,int n)
{
 psnterm	*t1,*t2;

 for ( t1=seq->terms+p1,t2=seq->terms+p2 ; n>0 ; n--,t1++,t2++ )
  {	if ( t1->type != t2->type || t1->major != t2->major )
		return(1);
  };	
 return(0);
}
static int psn_int_cpy(psn *seq,int p1,int p2,int n)
{
 memmove(seq->terms+p1,seq->terms+p2,n*sizeof(psnterm));
 return(0);
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

int psn_optimize(psn *pseq)
{
 psnterm 	*wseq;
 psn		*ncc;
 int		i,j,k,pn,ps,sp,stacklen,in,w; 
 istack		*stack;

 int	nopt,ofound;		/* number of optimized elements...	*/

 ncc=psn_alloc();
 if ( ncc==NULL )	return(psnerrno=PEALLOC);
 ncc->nseq=pseq->nseq;

 stacklen=PSNSTACKBLOCK;
 stack=(istack *)malloc(sizeof(istack)*stacklen);
 sp=0;
 stack[0].offset=0;
 stack[0].length=0;

 for ( wseq=pseq->terms ; wseq->type ; wseq++ )	/* already optimized!	*/
  {	if ( wseq->type==T_STACKVAR )
		return(psnerrno=0);
  } 

 nopt=0;
 
 for ( wseq=pseq->terms ; wseq->type ; wseq++ )
  { switch( wseq->type )
    { case T_CONST: case T_SCONST: case T_VAR:

	pn=stack[sp].offset;
	psn_append(ncc,wseq,1);
	pn++;
	stack[sp].length=1;
	sp++;
	stack[sp].offset=pn;
	stack[sp].length=0;

	break;

      case T_OP: case T_FN:

	in=wseq->minor;

	pn=stack[sp].offset;
	psn_append(ncc,wseq,1);

	sp-=in;	

	stack[sp].length++;
	for ( i=1 ; i<in ; i++ )
	 {	stack[sp].length+=stack[sp+i].length;	}
	pn=stack[sp].offset,
	ps=stack[sp].length;

	if ( ps>1 )
	 {	for ( i=0 ; i<nopt ; i++ )
		 {	if ( psn_int_cmp(ncc,stack[i].offset,pn,ps)==0 && 
			stack[i].length==ps )
				break;
		 }
		if ( i<nopt && nopt )	/* found... */
		 {	stack[sp].length=1;
			psn_remove(ncc,stack[sp].offset,ps);
			psn_append_term(ncc,T_STACKVAR,i,0,-1);
		 }
		else			/* not fount -> search! */
		 {	for ( j=stack[nopt].offset,ofound=0 ; 
			j+ps<=stack[sp].offset ; j++ )
			 {	if ( psn_int_cmp(ncc,j,stack[sp].offset,ps) != 0 )
					continue;

				if ( ! ofound )
				 {	w=stack[nopt].offset;
					psn_insert(ncc,stack[nopt].offset,ps);
					for ( k=nopt ; k<=sp ; k++ )
					 {	stack[k].offset+=ps;	}
					j+=ps;
					memmove(stack+nopt+1,stack+nopt,
						sizeof(istack)*(sp-nopt+1));
					sp++;
					if ( sp>=stacklen-2 )
					 {	stacklen+=PSNSTACKBLOCK;
						stack=(istack *)realloc(stack,
						sizeof(istack)*stacklen);
					 }
					stack[nopt].offset=w;
					stack[nopt].length=ps;
					psn_int_cpy(ncc,stack[nopt].offset,
					stack[sp].offset,ps);
					ofound=1,nopt++;
				 }
				psn_remove(ncc,j,ps-1);
				for ( k=0 ; k<=sp ; k++ )
				 {  if ( stack[k].offset>j )
					stack[k].offset-=ps-1;
				    if ( stack[k].offset<=j && j<stack[k+1].offset && k<sp )
					stack[k].length-=ps-1;
				 }
				ncc->terms[j].type=T_STACKVAR;
				ncc->terms[j].major=nopt-1;
				ncc->terms[j].minor=0;
				ncc->terms[j].cache=-1;
			 }
			if ( ofound )
			 {	stack[sp].length=1;
				psn_remove(ncc,stack[sp].offset,ps);
				psn_append_term(ncc,T_STACKVAR,nopt-1,0,-1);
			 }
		 }
	 }
			
	sp++;
	stack[sp].offset=stack[sp-1].offset+stack[sp-1].length;
	stack[sp].length=0;

	break;
     }
   if ( sp>=stacklen-2 )
    {	stacklen+=PSNSTACKBLOCK;
	stack=(istack *)realloc(stack,sizeof(istack)*stacklen);
    }

  }

 psn_swap_sequences(ncc,pseq);	
 psn_free(ncc);

 free(stack);

 return(psnerrno=0);
}

/*****************************************************************************/

static int psn_search_variable_code(int *varlist,int in,int major)
{
 if ( varlist==NULL )
  {	if ( major<in )
		return(major);
	else
		return(-1);
  }
 else
  {	int	i;
	for ( i=0 ; i<in ; i++ )
	 {	if ( varlist[i]==major )
			return(i);
	 }
	return(-1);
  }
}

psn * psn_replace(psn *pseq,int major,psn *replacement,int *varlist)
{
 psnterm 	*wseq,*wrep;
 psn		*ncc;
 int		i,pn,ps,rs,sp,stacklen,in,k;
 istack		*stack;
 double		d;

 ncc=psn_alloc();
 if ( ncc==NULL )
  {	psnerrno=PEALLOC;
	return(NULL);
  }
 ncc->nseq=pseq->nseq;

 if ( pseq->cons != NULL )
	psn_cons_append_more(ncc,pseq->cons,pseq->ncon);

 stacklen=PSNSTACKBLOCK;
 stack=(istack *)malloc(sizeof(istack)*stacklen);
 sp=0;
 stack[0].offset=0;
 stack[0].length=0;

 for ( wseq=pseq->terms ; wseq->type ; wseq++ )	/* already optimized */
  {	if ( wseq->type==T_STACKVAR )
	 {	psn_free(ncc);
		psnerrno=PEOPTIMIZED;
		return(NULL);
	 }
  } 

 for ( wseq=pseq->terms ; wseq->type ; wseq++ )
  { switch( wseq->type )
    { case T_CONST: case T_SCONST: case T_VAR:

	pn=stack[sp].offset;
	psn_append(ncc,wseq,1);
	pn++;
	stack[sp].length=1;
	sp++;
	stack[sp].offset=pn;
	stack[sp].length=0;

	break;

      case T_OP: case T_FN:

	in=wseq->minor;
	if ( sp<in )
	 {	psn_free(ncc);
		psnerrno=PEINVALID;
		return(NULL);
	 }
		
	if ( wseq->major == major )	/* replace by 'replacement'	*/
	 {	ps=0;
		for ( wrep=replacement->terms ; wrep->type ; wrep++ )
		 {	if ( wrep->type==T_VAR && (k=psn_search_variable_code(varlist,in,wrep->major)) >=0 )
			 {	k+=sp-in;
				psn_append_inseq(ncc,stack[k].offset,stack[k].length);
				ps+=stack[k].length;
			 }
			else if ( wrep->type != T_CONST )
			 {	psn_append(ncc,wrep,1);
				ps++;
			 }
			else if ( replacement->cons==NULL || wrep->major >= replacement->ncon )
			 {	free(stack);
				psn_free(ncc);
				psnerrno=PEINVALID;
				return(NULL);
			 }
			else
			 {	d=replacement->cons[wrep->major];
				for ( i=0 ; i<ncc->ncon ; i++ )
				 {	if ( ncc->cons[i]==d )	break;	}
				if ( i >= ncc->ncon )
				 {	i=ncc->ncon;
					psn_cons_append(ncc,d);
				 }
				psn_append_term(ncc,T_CONST,i,0,0);
				ps++;
			 }
		 }
		rs=0;
		for ( i=0 ; i<in ; i++ )
		 {	rs+=stack[sp-in+i].length;		}
		psn_remove(ncc,stack[sp-in].offset,rs);
		stack[sp-in].length=ps;
		sp-=in;
		sp++;
		stack[sp].offset=stack[sp-1].offset+stack[sp-1].length;
		stack[sp].length=0;
	 }
	else				/* otherwise do nothing		*/
	 {	psn_append(ncc,wseq,1);
		sp-=in;	
		stack[sp].length++;
		for ( i=1 ; i<in ; i++ )
		 {	stack[sp].length+=stack[sp+i].length;	}
		sp++;
		stack[sp].offset=stack[sp-1].offset+stack[sp-1].length;
		stack[sp].length=0;
	 }
 	break;
     }
   if ( sp>=stacklen-2 )
    {	stacklen+=PSNSTACKBLOCK;
	stack=(istack *)realloc(stack,sizeof(istack)*stacklen);
    }

  }

 free(stack);

 psnerrno=0;
 return(ncc);
}

/*****************************************************************************/

psn * psn_duplicate(psn *pseq)
{
 psn	 *ret;

 ret=psn_alloc();

 ret->nseq=pseq->nseq;

 if ( pseq->cons != NULL )
	psn_cons_append_more(ret,pseq->cons,pseq->ncon);

 psn_append(ret,pseq->terms,pseq->nterm);

 return(ret);
}

/*****************************************************************************/

int psn_concate(psn *pseq,psn *ps2)
{
 int	 i,type,major,minor,cache;
 psnterm *seq;
 
 for ( i=0 ; i<pseq->nterm ; i++ )
  {	if ( pseq->terms[i].type==T_STACKVAR )
		return(psnerrno=PEOPTIMIZED);
  }
 for ( i=0 ; i<ps2->nterm ; i++ )
  {	if ( ps2->terms[i].type==T_STACKVAR )	
		return(psnerrno=PEOPTIMIZED);
  }

 for ( seq=ps2->terms ; seq->type ; seq++ )
  {	type=seq->type;
	major=seq->major;
	minor=seq->minor;
	cache=seq->cache;
	if ( type==T_CONST )
	 {	for ( i=0 ; i<pseq->ncon ; i++ )
		 {	if ( pseq->cons[i]==ps2->cons[major] )	break;	}
		if ( i<pseq->ncon && pseq->ncon>0 )
		 {	major=i;		}
		else
		 {	psn_cons_append(pseq,ps2->cons[major]);
			major=pseq->ncon-1;
		 }
	 }
	psn_append_term(pseq,type,major,minor,cache);
  }

 pseq->nseq+=ps2->nseq;
 return(psnerrno=0);
}

/*****************************************************************************/

static psn *psn_subextract(psn *pseq,int **rpc,int *rnopt)
{
 psnterm *wseq;
 psn	*ncc;
 int	i,pn,sp,in,nopt,*pc,stacklen; 

 ncc=psn_alloc();
 if ( ncc==NULL )
  {	psnerrno=PEALLOC;
	return(NULL);
  }
 if ( pseq->cons != NULL )
	psn_cons_append_more(ncc,pseq->cons,pseq->ncon);

 pc=*rpc;
 stacklen=PSNSTACKBLOCK;
 sp=0;pc[0]=pc[1]=0;
 
 nopt=0;
 for ( wseq=pseq->terms ; wseq->type ; wseq++ )
  {	switch( wseq->type )
	 {   case T_CONST: case T_SCONST: case T_VAR:

		pn=pc[2*sp];
		psn_append(ncc,wseq,1);
		pn++;
		pc[2*sp+1]=1;

		sp++;
		pc[2*sp]=pn;
		pc[2*sp+1]=0;

		break;

	     case T_STACKVAR:

		i=wseq->major;
		psn_append_inseq(ncc,pc[2*i],pc[2*i+1]);
		if ( i+1>nopt )	nopt=i+1;

		pc[2*sp+1]=pc[2*i+1];
		
		sp++;
		pc[2*sp]=pc[2*(sp-1)]+pc[2*(sp-1)+1];
		pc[2*sp+1]=0;
		break;

	     case T_OP: case T_FN:

		in=wseq->minor;

		pn=pc[2*sp];
		psn_append(ncc,wseq,1);

		sp-=in;	

		pc[2*sp+1]++;
		for ( i=1 ; i<in ; i++ )
		 {	pc[2*sp+1]+=pc[2*(sp+i)+1];	}

		sp++;
		pc[2*sp]=pc[2*(sp-1)]+pc[2*(sp-1)+1];
		pc[2*sp+1]=0;

		break;
	 }
	if ( sp>=stacklen-2 )
	 {	stacklen+=PSNSTACKBLOCK;
		pc=(int *)realloc(pc,sizeof(int)*2*stacklen);
	 }
  }

 *rnopt=nopt;
 *rpc=pc;

 return(ncc);

}
psn *psn_extract(psn *pseq,int n)
{
 int	i,nopt;
 psn	*ncc;
 int	*pc;

 pc=(int *)malloc(sizeof(int)*2*PSNSTACKBLOCK);

 if ( n>=pseq->nseq )
  {	psnerrno=PEINVALID;
	return(NULL);
  }

 ncc=psn_subextract(pseq,&pc,&nopt);
 if ( ncc==NULL )
  {	/* psnerrno=PE... */ /*(already set by psn_subextract())*/
	return(NULL);
  }
 i=n+nopt;
 psn_remove(ncc,pc[2*(i+1)],ncc->nterm-pc[2*(i+1)]);
 psn_remove(ncc,0,pc[2*n]);
 ncc->nseq=1;

 free(pc);

 return(ncc);
}

psn *psn_full_extract(psn *pseq)
{
 int	nopt;
 psn	*ncc;
 int	*pc;

 pc=(int *)malloc(sizeof(int)*2*PSNSTACKBLOCK);

 ncc=psn_subextract(pseq,&pc,&nopt);
 if ( ncc==NULL )
  {	/* psnerrno=PE... */ /*(already set by psn_subextract())*/
	return(NULL);
  }

 psn_remove(ncc,0,pc[2*nopt]);
 ncc->nseq=pseq->nseq;

 free(pc);

 return(ncc);
}

/*****************************************************************************/

typedef struct
 {	char	*string;
	int	precedency;
	int	affix;
 } psnstringstack;

static int psn_convert_add_char(char **ret,int *rp,int *rsize,int ch)
{
 (*ret)[*rp]=ch;
 (*rp)++;
 if ( (*rp)>=(*rsize) )
	*ret=(char *)realloc(*ret,(*rsize)+BLOCKSIZE),
	*rsize+=BLOCKSIZE;
 return(0);
}
static int psn_convert_add_string(char **ret,int *rp,int *rsize,char *str)
{
 for ( ; *str ; str++ )
	psn_convert_add_char(ret,rp,rsize,*str);
 return(0); 
}

static char *psn_convert_replace_single(char *str,char *rep)
{
 char	*ret;
 int	rsize,rp;

 rsize=BLOCKSIZE,rp=0; 
 ret=(char *)malloc(rsize);
 while ( *str )
  {	if ( *str != '#' )
	 {	psn_convert_add_char(&ret,&rp,&rsize,*str);
		str++;
	 }
	else
	 {	str++;
		if ( *str != '0' )
		 {	psn_convert_add_char(&ret,&rp,&rsize,*str);
			str++;
		 }
		else
		 {	psn_convert_add_string(&ret,&rp,&rsize,rep);
			str++;
		 }
	 }
  };
 psn_convert_add_char(&ret,&rp,&rsize,0);
 return(ret);
}

static char *psn_convert_replace_more(char *str,psnstringstack *first,
char *parexpr,int precedency)
{
 char	*ret,*wexpr;
 int	rsize,rp,to_paren,to_paren_equ,n;
 int	is_term_first,is_term_last,affix;

 rsize=BLOCKSIZE,rp=0; 
 ret=(char *)malloc(rsize);

 is_term_first=1;

 while ( *str )
  {	if ( *str != '#' )
	 {	psn_convert_add_char(&ret,&rp,&rsize,*str);
		str++;
	 }
	else if ( *(str+1) == '#' )
	 {	psn_convert_add_char(&ret,&rp,&rsize,'#');
		str++;str++;
	 }
	else
	 {	str++;to_paren=0,to_paren_equ=0;
		if ( *str == '(' )	to_paren=1,str++;
		else if ( *str == '[' )	to_paren=1,to_paren_equ=1,str++;

		n=*str-'1',str++;
		if ( n<0 || n>8 )
		 {	free(ret);
			psnerrno=PEINVALID;
			return(NULL);
 		 }
		if ( to_paren )
		 {	if ( ( *str != ')' && ! to_paren_equ ) ||
			( *str != ']' && to_paren_equ ) )
			 {	free(ret);
				psnerrno=PEINVALID;
				return(NULL);
			 }
			str++;
		 }
	
		if ( *str==0 )	is_term_last=1;
		else		is_term_last=0;

		affix=first[n].affix;

		if ( to_paren && ( 
		( first[n].precedency <  precedency && ! to_paren_equ ) ||
		( first[n].precedency <= precedency && to_paren_equ ) ||
		( affix==TO_PREFIX && ! is_term_first ) ||
		( affix==TO_SUFFIX && ! is_term_last ) ) ) 
		 {	wexpr=psn_convert_replace_single(parexpr,
				first[n].string);
			psn_convert_add_string(&ret,&rp,&rsize,wexpr);
			free(wexpr);
		 }
		else
			psn_convert_add_string(&ret,&rp,&rsize,first[n].string);
	 }
	is_term_first=0;
  };
 psn_convert_add_char(&ret,&rp,&rsize,0);

 return(ret);
}

static char *psn_convert_to_const(double d,char *dformat)
{
 char	buff[128],*ret;
 if ( dformat==NULL )	sprintf(buff,"%g",d);
 else			sprintf(buff,dformat,d);
 ret=(char *)malloc(strlen(buff)+1);
 strcpy(ret,buff);
 return(ret);
}
static char *psn_convert_to_symbol(int major,psnsym **symtable)
{
 char	*name,*ret;
 name=psn_symtable_search(symtable,major,T_VAR,T_VAR);
 if ( name==NULL )	return(NULL);
 ret=(char *)malloc(strlen(name)+1);
 strcpy(ret,name);
 return(ret);
}

static int psn_convert_stack_cleanup(psnstringstack *chstack,int n)
{
 while ( n )
  {	free(chstack->string);
	n--,chstack++;
  };
 return(0);
}

char *psn_convert_symbolic(psn *seq,psnprop *plist,psnsymeval *syt,
			   psnsym **sym,char *dformat)
{
 psnterm	*wseq;
 psnprop	*wp;
 int	 	sp,argnum,i,major,stacklen;
 psnsymeval	*wsy,*sypar;
 char		*wstr,*ret;
 double		dd;
 psnstringstack	*chstack;

 if ( seq->nseq != 1 )
  {	psnerrno=PEMULTI;
	return(NULL);
  }

 for ( wseq=seq->terms ; wseq->type ; wseq++ )
  {	if ( wseq->type==T_STACKVAR )
	 {	psnerrno=PEOPTIMIZED;
		return(NULL);
	 }
  }

 for ( i=0,sypar=NULL ; ! ( syt[i].major==0 && syt[i].string==NULL ) ; i++ )
  {	if ( syt[i].major==0 )
	 {	sypar=&syt[i];
		break;
	 }
  }
 if ( sypar==NULL )
  {	psnerrno=PENOTFOUND;
	return(NULL);
  }

 stacklen=PSNSTACKBLOCK;
 chstack=(psnstringstack *)malloc(sizeof(psnstringstack)*stacklen);
 sp=0;

 for ( wseq=seq->terms ; wseq->type ; wseq++ )
  {	switch ( wseq->type )
	 {   case T_VAR:
		wstr=psn_convert_to_symbol(wseq->major,sym);
		if ( wstr==NULL )
		 {	psn_convert_stack_cleanup(chstack,sp);
			free(chstack);
			psnerrno=PENOTFOUND;
			return(NULL);
		 }
		chstack[sp].string=wstr;
		chstack[sp].precedency=PSN_MAX_PREC;
		chstack[sp].affix=0;
		sp++;
		break;
	     case T_SCONST: case T_CONST:
		if ( wseq->type==T_SCONST )	dd=(double)wseq->major;
		else				dd=seq->cons[wseq->major];

		wstr=psn_convert_to_const(dd,dformat);
		if ( dd<0.0 )
		 {	chstack[sp].string=
				psn_convert_replace_single(sypar->string,wstr);
			free(wstr);
		 }
		else
			chstack[sp].string=wstr;
			
		chstack[sp].precedency=PSN_MAX_PREC;
		chstack[sp].affix=0;

		sp++;

		break;
	     case T_OP: case T_FN:
		wp=psn_prop_search(plist,wseq->major);
		if ( wp==NULL )
		 {	psn_convert_stack_cleanup(chstack,sp);
			free(chstack);
			psnerrno=PENOTFOUND;
			return(NULL);
		 }
		argnum=wp->argnum;

		/* Currently, functions with arbitrary number of arguments   */
		/* (argnum<0) are not allowed here: the expected	     */
		/* number of arguments 'wp->argnum' should be equal to	     */
		/* the number of arguments stored in the auxiliary field     */
		/* 'wseq->minor'. However, if such desired number of	     */
		/* arguments cannot be found in the stack (thus, sp<argnum), */
		/* the conversion also should be stopped with an error.	     */

		if ( argnum > sp || wseq->minor != argnum ) 
		 {	psn_convert_stack_cleanup(chstack,sp);
			free(chstack);
			psnerrno=PEINVALID;
			return(NULL);
		 }

		for ( wsy=syt ; ! ( wsy->major == 0 && wsy->string == NULL ) ;
		 wsy++ )
		 {	if ( wsy->major == wseq->major )	break;	}
		if ( wsy->major == 0 && wsy->string == NULL )
		 {	psn_convert_stack_cleanup(chstack,sp);
			free(chstack);
			psnerrno=PENOTFOUND;
			return(NULL);
		 }

		major=wsy->major;
		wstr=psn_convert_replace_more(wsy->string,&chstack[sp-argnum],
		sypar->string,wp->precedency);
		if ( wstr==NULL )
		 {	psn_convert_stack_cleanup(chstack,sp);
			free(chstack);
			psnerrno=PENOTFOUND;
			return(NULL);
		 }
		for ( i=0 ; i<argnum ; i++ )
		 {	sp--;
			free(chstack[sp].string);
		 }

		chstack[sp].string=wstr;
		if ( wsy->strength )
			chstack[sp].precedency=PSN_MAX_PREC,
			chstack[sp].affix=0;
		else
			chstack[sp].precedency=wp->precedency;
			chstack[sp].affix=wsy->affixation;
		sp++;

		break;
	     default:
		psn_convert_stack_cleanup(chstack,sp);
		free(chstack);
		psnerrno=PEINVALID;
		return(NULL);
	 }
	if ( sp>=stacklen )
	 {	stacklen+=PSNSTACKBLOCK;
		chstack=(psnstringstack *)realloc(chstack,
		sizeof(psnstringstack)*stacklen);
	 }

  }

 if ( sp != 1 )
  {	psn_convert_stack_cleanup(chstack,sp);
	psnerrno=PEINVALID;
	free(chstack);
	return(NULL);
  }

 ret=chstack[0].string;
 free(chstack);
 return ( ret );
}	

/*****************************************************************************/

typedef struct
 {	int	op;
	int	last;
	int	prev;
	int	block;
	int	size;
 } psntree;

typedef struct
 {	int	major;
	short	*pre;
	int	prelen;
	short	*post;
	int	postlen;
	psntree	*tree;
	int	treepos;
 } psnsimpcache;

typedef struct 
 {	psnsimpcache	*psc;
	int		ns;
	psntree		*alltree;
 } psnsimpinfo;

typedef struct
 {	int	block;
	int	size;
 } psnsimpexp;

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

static int psn_build_simp_tree(short *srule,psntree *tree,psnprop *plist)
{
 short		*wseq;
 psntree	*wt;
 psnprop	*wp;
 int		*opstack,sp,stacklen;
 int		i,j,k,l,n,argnum,size;

 stacklen=PSNSTACKBLOCK;
 opstack=(int *)malloc(sizeof(int)*stacklen);
 sp=0;

 for ( wseq=srule,n=0,wt=tree ; *wseq ; wseq++,n++,wt++ )
  {	if ( *wseq < 0 )
	 {	opstack[sp]=n;
		sp++;
		wt->op=-1;
		wt->last=-1;
		wt->prev=-1;
		wt->block=n;
		wt->size=1;
	 }
	else
	 {	wp=psn_prop_search(plist,*wseq);
		if ( wp==NULL )	
		 {	free(opstack);
			return(1);
		 }
		argnum=wp->argnum;
		size=1;
		if ( argnum>0 )
		 {	k=wt->last=opstack[sp-1];
			for ( i=argnum-1,j=sp-1 ; i>0 ; i--,j-- )
			 {	l=tree[k].prev=opstack[j-1];
				tree[k].op=n;
				size+=tree[k].size;
				k=l;
			 }
			tree[k].prev=-1;
			tree[k].op=n;
			size+=tree[k].size;
			wt->block=tree[k].block;
		 }
		else
		 {	wt->block=n;
			wt->last=-1;
		 }
		wt->size=size;
		wt->op=-1;
		wt->prev=-1;
		sp-=argnum;
		opstack[sp]=n;
		sp++;
	 }
	if ( sp>=stacklen )
	 {	stacklen+=PSNSTACKBLOCK;
		opstack=(int *)realloc(opstack,sizeof(int)*stacklen);
	 }
  }

 free(opstack);
 return(0);
}

static int psn_build_simp_info(psnsimpinfo *psi,short **rules,psnprop *plist)
{
 int		i,j,k,ns;
 short		*rlist;
 psnsimpcache	*psc;
 psntree	*alltree;
 int		alltreelen;

 for ( ns=0 ; rules[ns] != NULL ; ns++ ) ;
 psc=(psnsimpcache *)malloc(sizeof(psnsimpcache)*ns);
 alltree=NULL;
 alltreelen=0;
 for ( i=0 ; i<ns ; i++ )
  {	rlist=rules[i];
	psc[i].pre=rlist;
	for ( j=0 ; *rlist ; rlist++,j++ )	;
	psc[i].prelen=j;
	psc[i].major=psc[i].pre[j-1];
	rlist++;
	psc[i].post=rlist;
	for ( j=0 ; *rlist ; rlist++,j++ )	;
	psc[i].postlen=j;

	psc[i].treepos=j=alltreelen;
	alltreelen+=psc[i].prelen;
	if ( alltree==NULL )
		alltree=(psntree *)malloc (sizeof(psntree)*alltreelen);
	else
		alltree=(psntree *)realloc(alltree,sizeof(psntree)*alltreelen);
	k=psn_build_simp_tree(rules[i],&alltree[j],plist);
	if ( k )
	 {	free(alltree);
		free(psc);
		return(1);
	 }	
  }

 for ( i=0 ; i<ns ; i++ )
  {	psc[i].tree=&alltree[psc[i].treepos];		}

 psi->ns=ns;
 psi->psc=psc;
 psi->alltree=alltree;
 return(0);
}

static int psn_drop_simp_info(psnsimpinfo *psi)
{
 if ( psi->psc != NULL )	free(psi->psc);
 if ( psi->alltree != NULL )	free(psi->alltree);
 return(0);
}

static int psn_build_tree(psn *seq,psntree *tree)
{
 psnterm	*wseq;
 psntree	*wt;
 int		*opstack,sp,stacklen;
 int		i,j,k,l,n,argnum,size;

 stacklen=PSNSTACKBLOCK;
 opstack=(int *)malloc(sizeof(int)*stacklen);
 sp=0;

 for ( wseq=seq->terms,n=0,wt=tree ; wseq->type ; wseq++,n++,wt++ )
  {	switch ( wseq->type )
	 {   case T_VAR: case T_CONST: case T_SCONST:
		opstack[sp]=n;
		sp++;
		wt->op=-1;
		wt->last=-1;
		wt->prev=-1;
		wt->block=n;
		wt->size=1;
		break;
	     case T_FN : case T_OP:
		argnum=wseq->minor;
		size=1;
		if ( argnum>0 )
		 {	k=wt->last=opstack[sp-1];
			for ( i=argnum-1,j=sp-1 ; i>0 ; i--,j-- )
			 {	l=tree[k].prev=opstack[j-1];
				tree[k].op=n;
				size+=tree[k].size;
				k=l;
			 }
			tree[k].prev=-1;
			tree[k].op=n;
			size+=tree[k].size;
			wt->block=tree[k].block;
		 }
		else
		 {	wt->block=n;
			wt->last=-1;
		 }
		wt->size=size;
		wt->op=-1;
		wt->prev=-1;
		sp-=argnum;
		opstack[sp]=n;
		sp++;
		break;
	 }

	if ( sp>=stacklen )
	 {	stacklen+=PSNSTACKBLOCK;
		opstack=(int *)realloc(opstack,sizeof(int)*stacklen);
	 }
  }

 free(opstack);
 return(0);
}

static int psn_simp_compare(psn *seq,psntree *tree,int t,
			    psnsimpexp *psex,
			    short *pre,int cp,psntree *stree)
{
 int		w,k,wt,ws;
 psnterm	*terms;

 terms=seq->terms;

 w=pre[cp];

 if ( w>0 )
  {	if ( ! ( terms[t].type == T_OP || terms[t].type == T_FN ) )
		return(1);
	if ( terms[t].major != w )
		return(1);
	wt=tree[t].last,
	ws=stree[cp].last;
	while ( ws>=0 )
	 {	if ( wt<0 )	return(1);
		k=psn_simp_compare(seq,tree,wt,psex,pre,ws,stree);	
		if ( k )	return(1);
		wt=tree [wt].prev;
		ws=stree[ws].prev;
	 };
	return(0);
  }
 else if ( w<=-2*SS_REF )
  {	int	major;

	major=terms[t].major;
	w=-2*SS_REF-w;

	if ( terms[t].type==T_SCONST && major==w )
		return(0);
	else if ( terms[t].type==T_CONST && seq->cons[major]==(double)w )
		return(0);
	else
		return(1);
  }
 else if ( w<0 )
  {	w=(-w)-1;
	if ( psex[w].block<0 )
	 {	psex[w].block=tree[t].block;
		psex[w].size =tree[t].size ;
		return(0);
	 }
	else
	 {	int	block,size,wblock,wsize;
		block=psex[w].block,
		size =psex[w].size;
		wblock=tree[t].block,
		wsize =tree[t].size;
		if ( size != wsize )	return(1);
		k=psn_int_cmp(seq,block,wblock,size);
		return(k);
	 }
  }
 else /* w==0, this should not happen. */
	return(-1);
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifdef __DEBUG
static void __psn_print_tree(psntree *tree,int n)
{
 int	i;
 FILE	*fw;
 fw=stderr;
 fprintf(fw,"op:   ");
 for ( i=0 ; i<n ; i++ ) { fprintf(fw,"%2d ",tree[i].op);   } fprintf(fw,"\n");
 fprintf(fw,"last: ");
 for ( i=0 ; i<n ; i++ ) { fprintf(fw,"%2d ",tree[i].last); } fprintf(fw,"\n");
 fprintf(fw,"prev: ");
 for ( i=0 ; i<n ; i++ ) { fprintf(fw,"%2d ",tree[i].prev); } fprintf(fw,"\n");
 fprintf(fw,"block:");
 for ( i=0 ; i<n ; i++ ) { fprintf(fw,"%2d ",tree[i].block);} fprintf(fw,"\n");
 fprintf(fw,"size: ");
 for ( i=0 ; i<n ; i++ ) { fprintf(fw,"%2d ",tree[i].size); } fprintf(fw,"\n");
}
#endif

int psn_simplify(psn *seq,short **rules,psnprop *plist,psnfunct *flist)
{
 psntree	*wtree;
 psnterm	*wterm;
 int		i,j,k,t;
 int		argnum,major,rfound,isize,osize,oblock,olen,w,wtreelen;
 psnsimpinfo	psi_data,*psi=&psi_data;
 psnsimpexp	psex[2*SS_REF];
 psnprop	*wp;
 short		*post,*ws;
 double		consstack[2*SS_REF];

 for ( wterm=seq->terms ; wterm->type ; wterm++ )
  {	if ( wterm->type==T_STACKVAR )
		return(psnerrno=PEOPTIMIZED);
  }

 k=psn_build_simp_info(psi,rules,plist);
 if ( k ) 	return(psnerrno=PEINVALID);

 wtreelen=seq->nterm;
 wtree=(psntree *)malloc(sizeof(psntree)*wtreelen);
 psn_build_tree(seq,wtree);

#ifdef __DEBUG
 __psn_print_tree(wtree,wtreelen);
#endif

 for ( t=0 ; seq->terms[t].type ; t++ )
  {	
	if ( ! ( seq->terms[t].type==T_OP || seq->terms[t].type==T_FN ) )
		continue;

	if ( seq->terms[t].minor>0 && flist != NULL )
	 {	int	is_all_const;
		double	*cnst,dret;
		short	sret;

		is_all_const=1;
		for ( j=wtree[t].last ; j>=0 && is_all_const; j=wtree[j].prev )
		 {	if ( ! ( wtree[j].size==1 && (
			seq->terms[j].type==T_SCONST || 
			seq->terms[j].type==T_CONST ) ) )
				is_all_const=0;
		 }
		if ( is_all_const )	/* All operands are constant... */
		 {	argnum=seq->terms[t].minor;
			cnst=consstack;
			t-=argnum;
			for ( i=0 ; i<argnum ; i++ )
			 {	major=seq->terms[t+i].major;
				if ( seq->terms[t+i].type==T_SCONST )
					cnst[i]=(double)major;
				else
					cnst[i]=seq->cons[major];
			 }
			i=psn_cache_search(seq->terms[t+argnum].major,flist);
			if ( i<0 )
			 {	free(wtree);
				psn_drop_simp_info(psi);
				return(psnerrno=PENOTFOUND);
			 }
			k=flist[i].funct(cnst+argnum);
			if ( k )
			 {	free(wtree);
				psn_drop_simp_info(psi);
				return(psnerrno=PENUMERICAL);
			 }
			dret=cnst[0];sret=(short)dret;
			psn_remove(seq,t,argnum);
			if ( dret != (double)sret )
			 {	major=psn_cons_append(seq,dret);
				seq->terms[t].type =T_CONST;
				seq->terms[t].major=major;
			 }
			else
			 {	seq->terms[t].type =T_SCONST;
				seq->terms[t].major=sret;
			 }
			seq->terms[t].minor=0;
			seq->terms[t].cache=0;
			psn_build_tree(seq,wtree);

			continue;
		 }
	 }

	major=seq->terms[t].major;
	for ( i=0,rfound=-1 ; i<psi->ns && rfound<0 ; i++ )
	 {	
		if ( major != psi->psc[i].major )
			continue;

		for ( j=0 ; j<2*SS_REF ; j++ )
		 {	psex[j].block=-1;
			psex[j].size =-1;
		 }


		if ( ! psn_simp_compare(seq,wtree,t,psex,
		psi->psc[i].pre,psi->psc[i].prelen-1,psi->psc[i].tree) )
			rfound=i;
	 }
	if ( rfound < 0 )
		continue;

	if ( psi->psc[rfound].postlen==0 )
	 {	free(wtree);
		psn_drop_simp_info(psi);
		return(psnerrno=PEANALYTICAL);
	 }

#ifdef	__DEBUG
	fprintf(stderr,"t=%d, rfound=%d\n",t,rfound);
#endif

	post=psi->psc[rfound].post;

	isize=0;
	osize =wtree[t].size;
	oblock=wtree[t].block;
	olen=seq->nterm;
	for ( ws=post ; *ws ; ws++ )
	 {	w=(*ws);
		if ( w>0 )
		 {	major=w;
			wp=psn_prop_search(plist,major);
			psn_append_term(seq,T_OP,major,wp->argnum,-1);
			isize++;
		 }
		else if ( w<=-(2*SS_REF) )
		 {	major=-w-(2*SS_REF);
			psn_append_term(seq,T_SCONST,major,0,0);
			isize++;
		 }
		else if ( w<0 )
		 {	j=-w-1;
			psn_append_inseq(seq,psex[j].block,psex[j].size);
			isize+=psex[j].size;
		 }
	 }

#ifdef	__DEBUG
	fprintf(stderr,	"oblock=%d, osize=%d, olen=%d, isize=%d\n",
			oblock,osize,olen,isize);
#endif

	if ( isize>osize )
	 {	psn_insert(seq,oblock,isize-osize);
		olen+=isize-osize;
		osize=isize;
	 }
	psn_int_cpy(seq,oblock,olen,isize);
	psn_remove(seq,olen,isize);
	if ( osize>isize )
	 {	psn_remove(seq,oblock+isize,osize-isize);	}
	if ( seq->nterm > wtreelen )
	 {	wtreelen=seq->nterm;
		wtree=(psntree *)realloc(wtree,sizeof(psntree)*wtreelen);
	 }
	psn_build_tree(seq,wtree);
	t=oblock-1;
  }

 free(wtree);
 psn_drop_simp_info(psi);

 psn_cache_clear(seq);

 return(psnerrno=0);
}

/*****************************************************************************/

int psn_register_function(psnsym **syms,psnprop **props,
	psnfunct **functs,psndiff **diffs,psnsymeval **symevals,
	int major,char *name,int argnum,int (*funct)(double *),
	short *diffrule,char *symevalstr)
{
 int		c;
 psnsym		*wsym;
 psnprop	*wprop;
 psnfunct	*wfunct;
 psndiff	*wdiff;
 psnsymeval	*wsymeval;

/*
 fprintf(stderr,"psn_register_function(): registering something with major=%d\n",major);
*/

 if ( props==NULL )
	return(1);

 if ( syms != NULL && name != NULL )
  {	if ( *syms==NULL )
	 	c=0;
	else
	 {	for ( c=0 ; (*syms)[c].type>0 ; )	c++;	}
	*syms=(psnsym *)realloc(*syms,sizeof(psnsym)*(c+2));
	wsym=(*syms)+c;
	wsym->type=T_FN;
	wsym->major=major;
	wsym->name=name;
	wsym->minor=argnum;
	memset((*syms)+c+1,0,sizeof(psnsym));
  }

 if ( *props==NULL )
 	c=0;
 else
   {	for ( c=0 ; (*props)[c].major>0 ; )	c++;	}
 *props=(psnprop *)realloc(*props,sizeof(psnprop)*(c+2));
 wprop=(*props)+c;
 wprop->major=major;
 wprop->argnum=argnum;
 wprop->precedency=0;
 wprop->associativity=0;
 memset((*props)+c+1,0,sizeof(psnprop));

 if ( functs != NULL && funct != NULL )
  {	if ( *functs==NULL )
	 	c=0;
	else
	  {	for ( c=0 ; (*functs)[c].major>0 ; )	c++;	}
	*functs=(psnfunct *)realloc(*functs,sizeof(psnfunct)*(c+2));
	wfunct=(*functs)+c;
	wfunct->major=major;
	wfunct->funct=funct;
	memset((*functs)+c+1,0,sizeof(psnfunct));
  }

 if ( diffs != NULL && diffrule != NULL )
  {	if ( *diffs==NULL )
	 	c=0;
	else
	  {	for ( c=0 ; (*diffs)[c].major>0 ; )	c++;	}
	*diffs=(psndiff *)realloc(*diffs,sizeof(psndiff)*(c+2));
	wdiff=(*diffs)+c;
	wdiff->major=major;
	wdiff->oplist=diffrule;
	memset((*diffs)+c+1,0,sizeof(psndiff));
  }

 if ( symevals != NULL && symevalstr != NULL )
  {	if ( *symevals==NULL )
	 	c=0;
	else
	  {	for ( c=0 ; (*symevals)[c].major>0 ; )	c++;	}
	*symevals=(psnsymeval *)realloc(*symevals,sizeof(psnsymeval)*(c+2));
	wsymeval=(*symevals)+c;
	wsymeval->major=major;
	wsymeval->string=symevalstr;
	wsymeval->strength=0;
	wsymeval->affixation=0;
	memset((*symevals)+c+1,0,sizeof(psnsymeval));
  }

/*
 fprintf(stderr,"psn_register_function(): registering something with major=%d\n",major);
*/

 return(0);
}

/*****************************************************************************/

psn *psn_diff_simplify(psn *seq,psnprop *plist,psndiff *dlist,
	short **slist,psnfunct *flist,int var)
{
 psn		*diff,*dup;
 int		ret;

 if ( slist != NULL )
  {	dup=psn_duplicate(seq);
	ret=psn_simplify(dup,slist,plist,flist);
	if ( ret || psn_test(dup) )
	 {	psn_free(dup);
		return(NULL);
	 }
	diff=psn_diff_int(dup,plist,dlist,var,1);
	psn_free(dup);
  }
 else
	diff=psn_diff_int(seq,plist,dlist,var,1);

 if ( diff==NULL )
	return(NULL);	/* psnerrno set by psn_diff();	*/

 if ( psn_test(diff) )
  {	psn_free(diff);	
	return(NULL);	/* psnerrno set by psn_diff();	*/
  }

 if ( slist != NULL )
  {	ret=psn_simplify(diff,slist,plist,flist);
	if ( ret || psn_test(diff) )
	 {	psn_free(diff);
		return(NULL);	/* psnerrno set by psn_diff();	*/
	 }
  }

 psnerrno=PEOK;
 return(diff);
}

/*****************************************************************************/

int * psn_argument_chain(psn *pseq)
{
 int		stack_automatic[PSNSTACKBLOCK],*stack_dynamic,*stack;
 int		i,j,nopt,narg,sp;
 psnterm	*term;
 int		*chain;
 int		stacklen,chainlen;

 stack=stack_automatic;
 stack_dynamic=NULL;
 stacklen=PSNSTACKBLOCK;
 chainlen=PSNSTACKBLOCK;

 chain=(int *)malloc(sizeof(int)*chainlen);

 nopt=0;sp=0;
 for ( term=pseq->terms,i=0 ; term->type ; term++,i++ )
  {	chain[i]=-1;
	switch( term->type )
	 {   case T_CONST: case T_VAR: case T_SCONST: case T_STACKVAR:
		stack[sp]=i,sp++;
	  	break;
	     case T_OP: case T_FN:
		narg=term->minor;
		if ( sp<narg )
		 {	free(chain);
			if ( stack_dynamic != NULL )
				free(stack_dynamic);
			psnerrno=PEINVALID;
			return(NULL);
		 }
		sp-=narg;
		for ( j=0 ; j<narg ; j++ )
		 	chain[stack[sp+j]]=i;
		stack[sp]=i,sp++;
		break;
	 }
	if ( i+1>=chainlen )
	 {	chainlen+=PSNSTACKBLOCK;
		chain=(int *)realloc(chain,sizeof(int)*chainlen);
	 }
	if ( sp>=stacklen )
	 {  if ( stacklen==PSNSTACKBLOCK )
	     {	int	i;
		stacklen+=PSNSTACKBLOCK;
		stack_dynamic=(int *)malloc(sizeof(int)*stacklen);
		for ( i=0 ; i<PSNSTACKBLOCK ; i++ )
		 {	stack_dynamic[i]=stack_automatic[i];	}
	     }
	    else
	     {	stacklen+=PSNSTACKBLOCK;
		stack_dynamic=(int *)realloc(stack_dynamic,
		sizeof(int)*stacklen);
	     }
	    stack=stack_dynamic;
	 }
  }

 if ( stack_dynamic != NULL )
	free(stack_dynamic);

 psnerrno=PEOK;
 return(chain);
}

/*****************************************************************************/
   
