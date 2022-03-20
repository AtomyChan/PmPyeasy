/*****************************************************************************/
/* spmatrix.c								     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* A small sparse matrix implementation. Written by Istvan Domsa (see 	     */
/* original copyright below), minor modifications (to be 'fi'-compatible)    */
/* were made by A. Pal (apal@szofi.elte.hu).				     */
/*****************************************************************************/

/* spmatrix.c : small sparse matrix implementation (domsa@konkoly.hu (2002)) */
/* written for the HAT Project */
/* spmatrix.c,v 5.5 2003/05/12 18:18:32 domsa Exp */

#include <stdlib.h>

#include "spmatrix.h"

/* this is a node in a sparse matrix */
typedef struct node
 {	int	row,col;	/* position in the matrix */
	int	value;		/* the integer value of the node */
	struct	node  * next_in_row;
	struct	node  * next_in_col;
 } node;

/* a 2D sparse matrix */
typedef struct
{	node	*elements;
	int    rows, cols;
	int    count;              /* how many elements in the matrix */
	sperr  err;
	node   *rwitness;
	node   *cwitness;
} sparse_matrix;

/* local functions to manipulate internal structures */
static void set_elem(node *this, int row, int col, int value,
		      node *nextrow, node *nextcol)
{
 if ( this != NULL )
  {	this->row   = row;
	this->col   = col;
	this->value = value;
	this->next_in_row = nextrow;
	this->next_in_col = nextcol;
  }
 return;
}

static node  *new_elem(int row,int col,int val,node *nextrow,node *nextcol)
{
 node	 *new=malloc(sizeof *new);

 if ( new != NULL )
	set_elem(new, row, col, val, nextrow, nextcol);
 return(new);
}

static node *find_node(sparse_matrix *sm,int row,int col)
{
 node	*header,*tmp,*nptr;

 if ( row > sm->rows || col > sm->cols )
	return(NULL);

 header = sm->elements;

 if ( col > row )
  {	nptr = row >= sm->rwitness->row ? sm->rwitness : header;
	while ( nptr->row < row )
		nptr = nptr->next_in_col;
	tmp = nptr->next_in_row;
	while ( tmp != nptr && tmp->col < col )
		tmp = tmp->next_in_row;
	if ( tmp->col == col )
		return(tmp);
  }
 else
  {	nptr = (col >= sm->cwitness->col ? sm->cwitness : header);
	while ( nptr->col < col )
		nptr = nptr->next_in_row;
	if ( !row )
		return(nptr);
	tmp = nptr->next_in_col;
	while ( tmp != nptr && tmp->row < row )
		tmp = tmp->next_in_col;
	if(tmp->row == row)
		return(tmp);
  }
 return(NULL);

}


/* create an empty sparse matrix */
spmatrix *spm_create(void)
{
 node		*head;
 sparse_matrix	*sm;
 spmatrix	*sm_cont=NULL;

 if (	(sm=malloc(sizeof *sm)) == NULL ||
	(sm_cont=malloc(sizeof *sm_cont)) == NULL ||
	(head=malloc(sizeof *head)) == NULL )
  {	free(sm);
	free(sm_cont);
	return(NULL);
  }

 set_elem(head,0,0,0,head,head);

 sm->elements=head;
 sm->rows    =0;
 sm->cols    =0;
 sm->count   =1;	/* head is the first element */
 sm->err     =SpNoErr;
 sm->rwitness=head;
 sm->cwitness=head;

 sm_cont->matrix=(void *) sm;

 return(sm_cont);
}

/* clean up a sparse matrix */
void spm_destroy(spmatrix *spm)
{
 sparse_matrix	*sm;
 node		*nd,*nd_next,*header,*first;
 int		count;

 if ( spm==NULL )	return;

 sm     =(sparse_matrix *) (spm->matrix);
 header =sm->elements;
 count  =sm->count;
 first  =header;
 nd_next=header->next_in_row;

 do
  {	nd = nd_next;
	if ( nd==first )
	 {	nd_next = first->next_in_col;
		if ( first != header )
			free(first);
		first   = nd_next;
		nd_next = first->next_in_row;
	 }
	else
	 {	nd_next = nd->next_in_row;
		free(nd);
	 }
  } while(--count);

 free(header);
 free(sm);
 free(spm);

 return;
}

/* get error status of a sparse matrix */
sperr spm_geterr(spmatrix *sm)
{
 if ( sm == NULL )	return(SpNull);
 else			return(((sparse_matrix *)(sm->matrix))->err);
}

/* get number of elements in the sparse matrix */
int  spm_numofelems(spmatrix  *sm)
{
 if ( sm == NULL )	return(0);
 else			return(((sparse_matrix *)(sm->matrix))->count);
}

/* set a node value. if it does not exist, create it */
void spm_setnode(spmatrix *spm,int row,int col,int value)
{
 node		*nd = NULL;
 sparse_matrix	*sm;

 if ( spm == NULL )	return;

 sm=(sparse_matrix *) spm->matrix;

 if ( row<0 || col<0 )
  {	sm->err = SpInvSub;
	return;
  }

 /* try to find the node */
 if ( col <= sm->cols && row <= sm->rows )
	nd = find_node(sm, row, col);

 /* create a new node and set its value */
 if ( nd == NULL )
  {	node  *prevrow, *prevcol, *tmp;
	node  *new = NULL, *colptr = NULL, *rowptr = NULL;

	if ( col > sm->cols )
	 {	int	ii,ccol;
		node	*befptr=find_node(sm, 0, sm->cols);

		for ( ii=sm->cols+1 ; ii<=col ; ii++ )
		 {	if ( (new=malloc(sizeof *new)) == NULL )
			 {	sm->err = SpMem;
				return;
			 }
			set_elem(new, 0, ii, 0, sm->elements, new);
			befptr->next_in_row = new;
			befptr = new;
			(sm->count)++;
		 }
		colptr   = new;
		sm->cols = col;

		/* set witness */
		ccol = col/2;
		while ( sm->cwitness->col<ccol )
			sm->cwitness = sm->cwitness->next_in_row;
	 }
	if ( row > sm->rows )
	 {	int	ii,rrow;
		node	*befptr=find_node(sm, sm->rows, 0);

		for ( ii=sm->rows+1 ; ii <= row ; ii++ )
		 {	if ( (new = malloc(sizeof *new)) == NULL )
			 {	sm->err = SpMem;
				return;
			 }
			set_elem(new, ii, 0, 0, new, sm->elements);
			befptr->next_in_col = new;
			befptr = new;
			(sm->count)++;
		 }
		rowptr   = new;
		sm->rows = row;

		/* set witness */
		rrow = row/2;
		while ( sm->rwitness->row<rrow )
			sm->rwitness = sm->rwitness->next_in_col;
	 }

	if ( ! row )		/* colptr already pointing to required node */
		(colptr->value)++;
	else if ( ! col )	/* rowptr already pointing to required node */
		(rowptr->value)++;
	else	/* find the column into which the new element will be inserted */
	 {	if ( colptr == NULL )
			colptr  = find_node(sm, 0, col);
		prevcol = colptr;
		tmp     = prevcol->next_in_col;
		while ( tmp != colptr && tmp->row<row )
		 {	prevcol = tmp;
			tmp     = tmp->next_in_col;
		 }
		/* find the row into which the new element will be inserted */
		if ( rowptr == NULL )
			rowptr  = find_node(sm, row, 0);
		prevrow = rowptr;
		tmp     = prevrow->next_in_row;
		while ( tmp != rowptr && tmp->col<col )
		 {	prevrow = tmp;
			tmp     = tmp->next_in_row;
		 }
		new=new_elem(row,col,1,prevrow->next_in_row, prevcol->next_in_col);
		if ( new == NULL )
		 {	sm->err = SpMem;
			return;
		 }
		prevcol->next_in_col = new;
		prevrow->next_in_row = new;
		(sm->count)++;
	 }
  }
 else	/* the node is already there: simply set its value */
	nd->value = value;
 return;
}

/* increase a node value by 1. if it does not exist, create it */
void spm_incrnode(spmatrix *spm,int row,int col)
{
 node		*nd=NULL;
 sparse_matrix	*sm;

 if ( spm == NULL )
	return;

 sm=(sparse_matrix *) spm->matrix;

 if ( row<0 || col<0 )
  {	sm->err = SpInvSub;
	return;
  }

 /* try to find the node */
 if ( col<=sm->cols && row<=sm->rows )
	nd = find_node(sm, row, col);

 /* create a new node and set its value */
 if ( nd == NULL )
  {	node  *prevrow, *prevcol, *tmp;
	node  *new = NULL, *colptr = NULL, *rowptr = NULL;

	if ( col > sm->cols )
	 {	int   ii, ccol;
		node  *befptr = find_node(sm, 0, sm->cols);

		for ( ii=sm->cols+1 ; ii<=col ; ii++ )
		 {	if ( (new = malloc(sizeof *new)) == NULL )
			 {	sm->err = SpMem;
				return;
			 }
			set_elem(new, 0, ii, 0, sm->elements, new);
			befptr->next_in_row = new;
			befptr = new;
			(sm->count)++;
		 }
		colptr   = new;
		sm->cols = col;

		/* set witness */
		ccol = col/2;
		while ( sm->cwitness->col<ccol )
			sm->cwitness = sm->cwitness->next_in_row;
	 }
	if ( row>sm->rows )
	 {	int	ii, rrow;
		node	*befptr = find_node(sm, sm->rows, 0);
		for ( ii=sm->rows+1 ; ii <= row ; ii++ )
		 {	if ( (new = malloc(sizeof *new)) == NULL )
			 {	sm->err = SpMem;
				return;
			 }
			set_elem(new, ii, 0, 0, new, sm->elements);
			befptr->next_in_col = new;
			befptr = new;
			(sm->count)++;
		 }
		rowptr   = new;
		sm->rows = row;

		/* set witness */
		rrow = row/2;
		while ( sm->rwitness->row<rrow )
			sm->rwitness = sm->rwitness->next_in_col;
	 }

	if ( ! row )		/* colptr already pointing to required node */
		(colptr->value)++;
	else if ( ! col )	/* rowptr already pointing to required node */
		(rowptr->value)++;
	else	/* find the column into which the new element will be inserted */
	 {	if ( colptr == NULL )
			colptr  = find_node(sm, 0, col);
		prevcol = colptr;
		tmp     = prevcol->next_in_col;
		while ( tmp != colptr && tmp->row < row )
		 {	prevcol = tmp;
			tmp     = tmp->next_in_col;
		 }
		/* find the row into which the new element will be inserted */
		if(rowptr == NULL)
		rowptr  = find_node(sm, row, 0);
		prevrow = rowptr;
		tmp     = prevrow->next_in_row;
		while ( tmp != rowptr && tmp->col < col )
		 {	prevrow = tmp;
			tmp     = tmp->next_in_row;
		 }
		new=new_elem(row,col,1,prevrow->next_in_row, prevcol->next_in_col);
		if ( new == NULL )
		 {	sm->err = SpMem;
			return;	
		 }
		prevcol->next_in_col = new;
		prevrow->next_in_row = new;
		(sm->count)++;
	 }
  }
 else	/* the node is already there: simply set its value */
	(nd->value)++;

 return;
}

/* increase a node value by the given value. if it does not exist, create it */
void  spm_addval(spmatrix *spm, int val, int row, int col)
{
 node		*nd = NULL;
 sparse_matrix	*sm;

 if ( spm == NULL )
	return;

 sm=(sparse_matrix *) spm->matrix;

 if ( row<0 || col < 0 )
  {	sm->err = SpInvSub;
	return;
  }

 /* try to find the node */
 if ( col <= sm->cols && row <= sm->rows )
	nd = find_node(sm, row, col);

 /* create a new node and set its value */
 if ( nd == NULL )
  {	node	*prevrow, *prevcol, *tmp;
	node	*new = NULL, *colptr = NULL, *rowptr = NULL;

	if ( col > sm->cols )
	 {	int	ii, ccol;
		node	*befptr=find_node(sm, 0, sm->cols);

		for ( ii=sm->cols+1 ; ii <= col ; ii++ )
		 {	if ( (new = malloc(sizeof *new)) == NULL )
			 {	sm->err = SpMem;
				return;
			 }
			set_elem(new, 0, ii, 0, sm->elements, new);
			befptr->next_in_row = new;
			befptr = new;
			(sm->count)++;
		 }
		colptr   = new;
		sm->cols = col;
	
		/* set witness */
		ccol = col/2;
		while ( sm->cwitness->col < ccol )
			sm->cwitness = sm->cwitness->next_in_row;
	 }
	if ( row > sm->rows )
	 {	int	ii,rrow;
		node  *befptr = find_node(sm, sm->rows, 0);
		for(ii = sm->rows + 1; ii <= row; ii++)
		 {	if ( (new = malloc(sizeof *new)) == NULL )
			 {	sm->err = SpMem;
				return;
			 }
			set_elem(new, ii, 0, 0, new, sm->elements);
			befptr->next_in_col = new;
			befptr = new;
			(sm->count)++;
		 }
		rowptr   = new;
		sm->rows = row;

		/* set witness */
		rrow = row/2;
		while ( sm->rwitness->row < rrow )
			sm->rwitness = sm->rwitness->next_in_col;
	 }

	if ( ! row )		/* colptr already pointing to required node */
		colptr->value += val;
	else if ( ! col )	/* rowptr already pointing to required node */
		rowptr->value += val;
	else	/* find the column into which the new element will be inserted */
	 {	if ( colptr == NULL )
			colptr  = find_node(sm, 0, col);
		prevcol = colptr;
		tmp     = prevcol->next_in_col;
		while ( tmp != colptr && tmp->row < row )
		 {	prevcol = tmp;
			tmp     = tmp->next_in_col;
		 }
		/* find the row into which the new element will be inserted */
		if ( rowptr == NULL )
		rowptr  = find_node(sm, row, 0);
		prevrow = rowptr;
		tmp     = prevrow->next_in_row;
		while ( tmp != rowptr && tmp->col < col )
		 {	prevrow = tmp;
			tmp     = tmp->next_in_row;
		 }
		new=new_elem(row,col,val,prevrow->next_in_row,prevcol->next_in_col);
		if ( new == NULL )
		 {	sm->err = SpMem;
			return;
		 }
		prevcol->next_in_col = new;
		prevrow->next_in_row = new;
		(sm->count)++;
	 }
  }
 
 else	/* the node is already there: simply set its value */
	nd->value += val;

 return;
}

/* return the value of a node */
int spm_getval(spmatrix *sm,int row,int col)
{
 node	*nd;

 if ( sm == NULL || row < 0 || col < 0 )
	return(0);

 nd=find_node((sparse_matrix *) (sm->matrix),row,col);

 return(nd==NULL ? 0:nd->value);
}

/* return a node or NULL */
spnode spm_getnode(spmatrix *spm, int row, int col)
{
 if ( spm == NULL || row < 0 || col < 0 )
	return(NULL);
 else	
	return((spnode)find_node((sparse_matrix *)(spm->matrix),row,col));
} 

/* get next non-zero value from the sparse-matrix */
spnode spm_getnextnode(spmatrix *spm, spnode snd, spnode *witness)
{
 node		*header, *nd, *first;
 sparse_matrix	*sm;

 if ( spm == NULL )
	return(NULL);

 sm    =(sparse_matrix *)(spm->matrix);
 header=sm->elements;

 if ( snd == NULL )
  {	nd       = header;
	first    = header;
	*witness = (spnode) header;
	if ( header->value != 0 )
		return((spnode) header);
  }
 else
  {	nd    = (node *) snd;
	first = (node *) *witness;
  }

 do
  {	nd = nd->next_in_row;
	if ( nd==first )
	 {	first   = first->next_in_col;

		/* no more elements left */
		if ( first == header )
			return  NULL;
		nd       = first;
		*witness = (spnode) first;
	 }
	if(nd->value)
	      return((spnode)nd);
  } while(nd != header);

 return(NULL);
}

/* returns column number of a node */
int spm_getcol(spnode snd)
{
 if ( snd == NULL )	return(-1);
 else			return(((node *) snd)->col);
}

/* returns row number of a node */
int spm_getrow(spnode snd)
{
 if ( snd == NULL )	return(-1);
 else			return(((node *) snd)->row);
}

/* returns value of a node */
int spm_getnodeval(spnode snd)
{
 if ( snd == NULL )	return(0);
 else			return(((node *) snd)->value);
}

/* set a node's value */
void spm_setnodeval(spnode snd, int value)
{
 if ( snd == NULL )	return;

 ((node *) snd)->value = value;
 return;
}
                   
