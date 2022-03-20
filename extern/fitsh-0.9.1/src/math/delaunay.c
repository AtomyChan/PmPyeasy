/*****************************************************************************/
/* delaunay.c								     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Functions related to Delaunay-triangulation. Written by Istvan Domsa (see */
/* original copyright below), minor modifications (to be 'fi'-compatible)    */
/* and extensions (see `make_neighbour_list()` and realted functions) 	     */
/* were made by A. Pal (apal@szofi.elte.hu).				     */
/*****************************************************************************/

/* triangulation.c : triangulate data points (domsa@konkoly.hu (2002)) */
/* written for the HAT Project */
/* triangulation.c,v 5.5 2003/05/13 16:09:44 domsa Exp */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "tpoint.h"
#include "delaunay.h"

/* some shortcuts */

#define SWAP(a,b,type) do { type tmp;tmp=a;a=b;b=tmp; } while(0) 

#define	VECTOR(p1,p2,u,v)	( (u)=(p2)->x-(p1)->x,(v)=(p2)->y-(p1)->y )

#define	CROSS_PROD_3P(p1,p2,p3)	(((p2)->x-(p1)->x)*((p3)->y-(p1)->y)-((p2)->y-(p1)->y)*((p3)->x-(p1)->x))

#define CROSS_PROD_2V(u1,v1,u2,v2)	( (u1)*(v2)-(v1)*(u2) )
#define DOT_PROD_2V(u1,v1,u2,v2)	( (u1)*(u2)+(v1)*(v2) )

#define	OTHER_POINT(e,p)	( (e)->orig==p ? (e)->dest  : (e)->orig  )

#define	PREV(e,p)		( (e)->orig==p ? (e)->oprev : (e)->dprev )
#define	NEXT(e,p)		( (e)->orig==p ? (e)->onext : (e)->dnext )

#define IDENT_REFS(e1,e2)	( (e1)==(e2) )

typedef enum
 {	right,
	left
 } side;

typedef struct edge    edge;
typedef struct dpoint  dpoint;

struct  edge
 {	dpoint	*orig;
	dpoint	*dest;
	edge	*onext,*oprev;
	edge	*dnext,*dprev;
	edge	*next,*prev;
 };

struct  dpoint
 {	int	ident;
	double  x,y;
	edge	*entry_pt;
 };

#define		LISTSIZE	128

typedef struct
 {	edge	*first;
	edge	*last;
 } tridata;

/* cross product of two points */
static void crossdot_prod(dpoint *p1,dpoint *p2,dpoint *p3,double *cross,double *dot)
{
 double	dx1,dx2,dy1,dy2;

 VECTOR(p1,p2,dx1,dy1);
 VECTOR(p1,p3,dx2,dy2);

 *cross=CROSS_PROD_2V(dx1,dy1,dx2,dy2);
 *dot  =  DOT_PROD_2V(dx1,dy1,dx2,dy2);

 return;
}

/* sort the vertices of the triangles */
void vertice_sort(double *sides,tpoint **vertices)
{
 if ( sides[0]<sides[1] )
  {	SWAP(sides[0],sides[1],double);
	SWAP(vertices[0],vertices[2],tpoint *);
  }
 if ( sides[1]<sides[2] )
  {	SWAP(sides[1],sides[2],double);
	SWAP(vertices[0],vertices[1],tpoint *);
	if ( sides[0]<sides[1] )
	 {	SWAP(sides[0],sides[1],double);
		SWAP(vertices[0],vertices[2],tpoint *);
	 }
  }  
 return;
}

static int tri_init(tridata *td)
{
 td->first=NULL;
 td->last=NULL;
 return(0);
}
static int tri_fini(tridata *td)
{
 edge	*e,*n;
 for ( e=td->first ; e != NULL ; )
  {	n=e->next;
	free(e);
	e=n;
  };
 return(0);
}
static edge * tri_edge_alloc(tridata *td)
{
 edge	*e;

 e=(edge *)malloc(sizeof(edge));
 if ( td->first==NULL )	
  {	td->first=e;		}
 
 e->prev=td->last;
 e->next=NULL;
 if ( e->prev != NULL )	e->prev->next=e;
 td->last=e;

 return(e);
}

static int tri_edge_free(tridata *td,edge *e)
{
 if ( e->prev != NULL )	e->prev->next=e->next;
 if ( e->next != NULL )	e->next->prev=e->prev;
 if ( td->first==e )	td->first=e->next;
 if ( td->last==e )	td->last=e->prev;
 free(e);
 return(0);
}

/* create an edge */
static edge *makeedge(tridata *td,dpoint *p1,dpoint *p2)
{
 edge	*ne;

 if ( p1==NULL || p2==NULL )
	return(NULL);

 ne=tri_edge_alloc(td);
/* ne=(edge *)malloc(sizeof(edge)); */

 ne->orig=p1;
 ne->dest=p2;

 ne->onext=ne;
 ne->oprev=ne;
 ne->dnext=ne;
 ne->dprev=ne;

 if ( p1->entry_pt==NULL )	p1->entry_pt=ne;
 if ( p2->entry_pt==NULL )	p2->entry_pt=ne;

 return(ne);
}

/* delete an edge */
static void deledge(tridata *td,edge *ed)
{
 dpoint  *p1,*p2;

 if ( ed==NULL )
	return;

 p1=ed->orig;
 p2=ed->dest;

 /* adjust entry points */
 if ( p1->entry_pt==ed )	p1->entry_pt=ed->onext;
 if ( p2->entry_pt==ed )	p2->entry_pt=ed->dnext;

 /* four edge links to change */
 if ( (ed->onext)->orig==p1 )	(ed->onext)->oprev=ed->oprev;
 else				(ed->onext)->dprev=ed->oprev;

 if ( (ed->oprev)->orig==p1 )	(ed->oprev)->onext=ed->onext;
 else				(ed->oprev)->dnext=ed->onext;

 if ( (ed->dnext)->orig==p2 )	(ed->dnext)->oprev=ed->dprev;
 else				(ed->dnext)->dprev=ed->dprev;

 if ( (ed->dprev)->orig==p2 )	(ed->dprev)->onext=ed->dnext;
 else				(ed->dprev)->dnext=ed->dnext;

 tri_edge_free(td,ed); 
/* free(ed); */
 return;
}

/* add an edge to a ring of edges */
static void spliceedge(tridata *td,edge *e1,edge *e2,dpoint *pp)
{
 edge	*next;

 /* e2 must be the unnattached edge and e1 the previous ccw edge to e2 */

 if ( e1->orig==pp )		next=e1->onext,e1->onext=e2;
 else				next=e1->dnext,e1->dnext=e2;

 if ( next->orig==pp )		next->oprev=e2;
 else				next->dprev=e2;

 if ( e2->orig==pp )		e2->onext=next,e2->oprev=e1;
 else				e2->dnext=next,e2->dprev=e1;

 return;
}

/* creates a new edge and adds it to two rings of edges */
/* p1 and p2 are the two vertices which are being joined */
/* e1 and e2 are the two edges associated with u and v respectively */
static edge *joinedge(tridata *td,edge *e1,dpoint *p1,edge *e2,dpoint *p2,side sd)
{
 edge  *enew;

 if ( (enew=makeedge(td,p1,p2))==NULL )
	return(NULL);

 if ( sd==left )
  {	if ( e1->orig==p1 )	spliceedge(td,e1->oprev,enew,p1);
	else			spliceedge(td,e1->dprev,enew,p1);
	spliceedge(td,e2,enew,p2);
  }
 else
  {	spliceedge(td,e1,enew,p1);
	if ( e2->orig==p2 )	spliceedge(td,e2->oprev,enew,p2);
	else			spliceedge(td,e2->dprev,enew,p2);
  }

 return(enew);
}

/* lower tangent of two triangulation */
static void  lower_tangent(tridata *td,edge *r_cw_l,dpoint *p1,edge *l_ccw_r,dpoint *p2,
				edge **l_lower,dpoint **org_l_lower,
				edge **r_lower,dpoint **org_r_lower)
{
 edge	*eleft,*eright;
 dpoint	*o_l,*o_r,*d_l,*d_r;

 eleft =r_cw_l;
 eright=l_ccw_r;
 o_l   =p1;
 d_l   =OTHER_POINT(eleft,p1);
 o_r   =p2;
 d_r   =OTHER_POINT(eright,p2);

 for ( ; ; )
  {	if ( CROSS_PROD_3P(o_l,d_l,o_r)>0.0 )
	 {	eleft=PREV(eleft,d_l);
		o_l  =d_l;
		d_l  =OTHER_POINT(eleft,o_l);
	 }
	else if ( CROSS_PROD_3P(o_r,d_r,o_l)<0.0 )
	 {	eright=NEXT(eright,d_r);
		o_r   =d_r;
		d_r   =OTHER_POINT(eright,o_r);
	 }
	else
		break;
  }

 *l_lower    =eleft;
 *r_lower    =eright;
 *org_l_lower=o_l;
 *org_r_lower=o_r;

 return;
}

/* merge two adjacent Delaunay triangulations into a single one */
static void mergehulls(tridata *td,edge *r_cw_l,dpoint *p1,edge *l_ccw_r,dpoint *p2,edge **l_tangent)
{
 double	c_p_l_cand,c_p_r_cand,d_p_l_cand,d_p_r_cand;
 double	cot_l_cand=0.0,cot_r_cand=0.0;
 int	above_l_cand,above_r_cand;
 edge	*l_lower,*r_lower,*base,*l_cand,*r_cand;
 dpoint	*org_r_lower,*org_l_lower,*org_base,*dest_base;
 dpoint	*dest_l_cand,*dest_r_cand;

 /* create first cross edge by joining lower common tangent */
 lower_tangent(td,r_cw_l,p1,l_ccw_r,p2,&l_lower,&org_l_lower,&r_lower,&org_r_lower);

 base     =joinedge(td,l_lower,org_l_lower,r_lower,org_r_lower,right);
 org_base =org_l_lower;
 dest_base=org_r_lower;

 /* we have to return lower tangent */
 *l_tangent=base;

 /* main merge loop */
 for ( ; ; )
  {	l_cand=NEXT(base,org_base);
	r_cand=PREV(base,dest_base);
	dest_l_cand=OTHER_POINT(l_cand,org_base);
	dest_r_cand=OTHER_POINT(r_cand,dest_base);

	crossdot_prod(dest_l_cand,org_base,dest_base,&c_p_l_cand,&d_p_l_cand);
	crossdot_prod(dest_r_cand,org_base,dest_base,&c_p_r_cand,&d_p_r_cand);

	above_l_cand=c_p_l_cand > 0.0;
	above_r_cand=c_p_r_cand > 0.0;

	if ( !above_l_cand && !above_r_cand )
		break;

	/* advance l_cand counter-clockwise, until the in_circle test fails */
	if ( above_l_cand )
	 {	double	c_p_next,d_p_next,cot_next;
		edge	*next;
		dpoint	*dest_next;

		cot_l_cand=d_p_l_cand/c_p_l_cand;

		for ( ; ; )
		 {	next     =NEXT(l_cand,org_base);
			dest_next=OTHER_POINT(next,org_base);

			crossdot_prod(dest_next,org_base,dest_base,&c_p_next,&d_p_next);

			if ( c_p_next <= 0.0 )		break;
			cot_next=d_p_next/c_p_next;
			if ( cot_next > cot_l_cand )	break;

			deledge(td,l_cand);	
			l_cand    =next;
			cot_l_cand=cot_next;
		 }
	 }
	/* advance r_cand clockwise, until the in_circle test fails */
	if ( above_r_cand )
	 {	double	c_p_prev,d_p_prev,cot_prev;
		edge    *prev;
		dpoint  *dest_prev;

		cot_r_cand=d_p_r_cand/c_p_r_cand;

		for ( ; ; )
		 {	prev     =PREV(r_cand,dest_base);
			dest_prev=OTHER_POINT(prev,dest_base);

			crossdot_prod(dest_prev,org_base,dest_base,&c_p_prev,&d_p_prev);

			if ( c_p_prev <= 0.0 )		break;
			cot_prev=d_p_prev/c_p_prev;
			if ( cot_prev > cot_r_cand )	break;

			deledge(td,r_cand);
			r_cand    =prev;
			cot_r_cand=cot_prev;
		 }
	 }
	/* add a cross edge from base to either l_cand or r_cand */
	dest_l_cand=OTHER_POINT(l_cand,org_base);
	dest_r_cand=OTHER_POINT(r_cand,dest_base);

	if ( !above_l_cand || (above_l_cand && above_r_cand && cot_r_cand<cot_l_cand) )
	 {	/* connect to the right */
		base     =joinedge(td,base,org_base,r_cand,dest_r_cand,right);
		dest_base=dest_r_cand;
	 }
	else
	 {	/* connect to the left */
		base    =joinedge(td,l_cand,dest_l_cand,base,dest_base,right);
		org_base=dest_l_cand; 
	 }
 }

 return;
}

/* Delaunay triangulation by divide and conquer */
static void divide(tridata *td,dpoint *points,int ileft,int iright,edge **l_ccw,edge **r_cw)
{
 int	cardinal,split;

 cardinal=iright-ileft+1;

 if ( cardinal==2 )		/* bottom of the recursion: make an edge */
  {	*l_ccw=makeedge(td,points+ileft,points+iright);
	*r_cw =*l_ccw;
  }
 else if ( cardinal==3 )	/* bottom of the recursion: make an edge or a triangle */
  {	double	cp;
	edge	*ea,*eb,*ec;
	dpoint	*pa=points+ileft;
	dpoint	*pb=pa+1;
	dpoint	*pc=points+iright;

	ea=makeedge(td,pa,pb);
	eb=makeedge(td,pb,pc);
	spliceedge(td,ea,eb,pb);

	cp=CROSS_PROD_3P(pa,pb,pc);

	if ( cp>0.0 )		/* make a triangle */
	 {	ec=joinedge(td,ea,pa,eb,pc,right);
		*l_ccw=ea;
		*r_cw =eb;
	 }
	else if ( cp<0.0 )	/* make a triangle */
	 {	ec=joinedge(td,ea,pa,eb,pc,left);
		*l_ccw=ec;
		*r_cw =ec;
	 }
	else			/* points are collinear */
	 {	*l_ccw=ea;
		*r_cw =eb;
	 }
  }
 else if ( cardinal > 3 )
  {	edge  *l_ccw_l,*r_cw_l,*l_ccw_r,*r_cw_r,*l_tangent;

	split=(ileft+iright)/2;
	divide(td,points,ileft,split,&l_ccw_l,&r_cw_l);
	divide(td,points,split+1,iright,&l_ccw_r,&r_cw_r);
	mergehulls(td,r_cw_l,points+split,l_ccw_r,points+split+1,&l_tangent);
	if ( l_tangent->orig==points+ ileft )	l_ccw_l=l_tangent;
	if ( l_tangent->dest==points+iright )	r_cw_r =l_tangent;
	*l_ccw=l_ccw_l;
	*r_cw =r_cw_r;
  }
 return;
}

/* sort points along the x-axis (for qsort()) */
static int sort_x(const void *p1,const void *p2)
{
 double	x1=((dpoint *) p1)->x;
 double	x2=((dpoint *) p2)->x;
 if ( x1==x2 )	return ( ((dpoint *)p2)->y<((dpoint *)p1)->y ? -1 : 1 );
 else		return ( x2<x1 ? 1 : -1 );
} 

/* turn point list to dpoints */
static dpoint * point2dpoint(int len,tpoint *points)
{
 int	ii;
 dpoint	*dpoints;

 if ( (dpoints=malloc(len*sizeof(*dpoints)) )==NULL )
	return(NULL);

 for( ii=0 ; ii<len ; ii++ )
  {	dpoints[ii].ident=ii;
	dpoints[ii].x=points[ii].xcoord;
	dpoints[ii].y=points[ii].ycoord;
	dpoints[ii].entry_pt=NULL;
  }
 return(dpoints);
}

static triangle *dp2tri(int len,tpoint *points,dpoint *dpoints,int *ntri)
{
 triangle	*triangles=NULL;
 edge		*e_start,*ep,*next;
 dpoint		*p1,*p2,*p3;
 int		ii,numtri;
 tpoint		*pn1,*pn2,*pn3;

 numtri=0;

 for ( ii=0 ; ii<len ; ii++ )
  {	p1=dpoints+ii;
	e_start=p1->entry_pt;
	ep=e_start;
	if ( e_start==NULL )		/* internal error */
	 {	free(triangles);
		return(NULL);
	 }
	do
	 {	p2=OTHER_POINT(ep,p1);

		if ( p1<p2 )
		 {	next=NEXT(ep,p1);
			p3=OTHER_POINT(next,p1);
			if ( p1<p3 )
			 {	if ( IDENT_REFS(NEXT(next,p3),PREV(ep,p2)) )
				 {	if ( p2>p3 )	SWAP(p2,p3,dpoint *);
					pn1=points+p1->ident;
					pn2=points+p2->ident;
					pn3=points+p3->ident;
					triangles=realloc(triangles,(numtri+1)*sizeof(*triangles));
					triangles[numtri].vertices[0]=pn1;
					triangles[numtri].vertices[1]=pn2;
					triangles[numtri].vertices[2]=pn3;
					triangles[numtri].trix=0.0;
					triangles[numtri].triy=0.0;
					numtri++;
				 }
			 }
	 	 }
		ep=NEXT(ep,p1);	/* next edge around p1 */
	 } while ( ! IDENT_REFS(ep,e_start) );
  }
 if ( ntri != NULL )	*ntri=numtri;
 return(triangles);
}

static int add_neig_point(tpoint ***rneigpoints,int *apnt,int *alen,tpoint *p)
{
 if ( *apnt >= *alen )
  {	(*rneigpoints)=(tpoint **)realloc((*rneigpoints),sizeof(tpoint *)*((*alen)+LISTSIZE));
	(*alen)+=LISTSIZE;
  }
 (*rneigpoints)[(*apnt)]=p;
 (*apnt)++;
 return(0);
}

static int make_neighbour_list(dpoint *dpoints,tpoint *tpoints,int n,int *neigindx,triinfo *ti)
{
 int	alen,apnt,i;
 dpoint	*d,*nd,*fd;
 tpoint	*t;
 edge	*e;

 apnt=alen=0;
 ti->neigpoints=NULL;
 for ( i=0,d=dpoints ; i<n ; i++,d++ )
  {	e=d->entry_pt;
	nd=fd=OTHER_POINT(e,d);
	neigindx[d->ident]=apnt;
	do
	 {	t=tpoints+nd->ident;
		add_neig_point(&ti->neigpoints,&apnt,&alen,t);
		e=NEXT(e,d);
		nd=OTHER_POINT(e,d);
	 } while ( nd != fd );
	add_neig_point(&ti->neigpoints,&apnt,&alen,NULL);
  };
 return(0);
}

int delaunay_triangulation(tpoint *points,int len,triangle **rtris,int *rntri,triinfo *ti)
{
 dpoint		*dpoints;
 triangle	*tris;
 int		ntri;
 edge		*l_cw,*r_ccw;
 tridata	td;

 if ( len < 3 || points==NULL )				return(1);
 if ( (dpoints=point2dpoint(len,points))==NULL )	return(1);

 qsort((void *)dpoints,len,sizeof(dpoint),sort_x);
 tri_init(&td);
 divide(&td,dpoints,0,len-1,&l_cw,&r_ccw);
 ntri=0;
 tris=dp2tri(len,points,dpoints,&ntri);

 if ( ti != NULL )
  {	int	i,*neigindx;
	neigindx=(int *)malloc(sizeof(int)*len);
	make_neighbour_list(dpoints,points,len,neigindx,ti);
	ti->neigs=(tpoint ***)malloc(sizeof(tpoint **)*len);
	for ( i=0 ; i<len ; i++ )
	 {	ti->neigs[i]=&ti->neigpoints[neigindx[i]];		}
	free(neigindx);	
  }

 tri_fini(&td); 
 free(dpoints);

 if ( rtris != NULL )	*rtris=tris;
 if ( rntri != NULL )	*rntri=ntri;

 return(0);
}

/*****************************************************************************/

int free_triangulation_info(triinfo *ti)
{
 if ( ti==NULL )
	return(-1);
 if ( ti->neigpoints != NULL )
  {	free(ti->neigpoints);
	ti->neigpoints=NULL;
  }
 if ( ti->neigs != NULL )
  {	free(ti->neigs);
	ti->neigs=NULL;
  }
 return(0);
}

/*****************************************************************************/

static int compare_triangle_vertices(const void *vp1,const void *vp2)
{
 triangle	*t1=(triangle *)vp1;
 triangle	*t2=(triangle *)vp2;
 return ( memcmp(&t1->vertices[0],&t2->vertices[0],3*sizeof(tpoint*)) );
}

#define		STATICNNEIG		16
#define		LISTBLOCK		32

int expand_triangulation(triinfo *ti,triangle **rtris,int *rntri,tpoint *points,int npoint,int level)
{
 triangle	*tris,*wtris;
 int		ntri,nwtri,atri,i,j,k,l,m,p,a,b,c;
 int		*ngs,nng,ang,nng0,nng1;
 tpoint		**tn;

 if ( ti==NULL || points==NULL )	return(1);

 atri=LISTBLOCK;
 tris=(triangle *)malloc(sizeof(triangle)*atri); 
 ntri=0;

 ang=LISTBLOCK;
 ngs=(int *)malloc(sizeof(int)*ang);

 for ( i=0 ; i<npoint ; i++ )
  {	tn=ti->neigs[i];
	nng=0;

	for ( j=0 ; tn[j] != NULL ; j++ )
	 {	k=tn[j]-points;
		/*if ( k<=i )	continue;*/
		if ( nng>=ang )
		 {	ang+=LISTBLOCK;
			ngs=(int *)realloc(ngs,sizeof(int)*ang);
		 }
		ngs[nng]=k;
		nng++;
	 }

	nng0=0;
	for ( l=1 ; l<level ; l++ )
	 {	nng1=nng;
		for ( m=nng0 ; m<nng1 ; m++ )
		 {	tn=ti->neigs[ngs[m]];
			for ( j=0 ; tn[j] != NULL ; j++ )
			 {	k=tn[j]-points;
				if ( k<=i )	continue;
				for ( p=0 ; p<nng1 ; p++ )
				 {	if ( k<=ngs[p] )	break;	}
				if ( p<nng1 )	continue;
				if ( nng>=ang )
				 {	ang+=LISTBLOCK;
					ngs=(int *)realloc(ngs,sizeof(int)*ang);
				 }
				ngs[nng]=k;
				nng++;	
			 }
		 }
		nng0=nng1;
	 }

	for ( k=0 ; k<nng ; k++ )
	 {	/*if ( ngs[k]>=i )	continue;*/
		for ( l=0 ; l<nng ; l++ )
		 {	if ( ngs[l]<=ngs[k] )	continue;
			if ( ntri>=atri )
			 {	atri+=LISTBLOCK;
				tris=(triangle *)realloc(tris,sizeof(triangle)*atri);
			 }
			a=i,
			b=ngs[k];
			c=ngs[l];
			if ( a>b )	p=a,a=b,b=p;
			if ( b>c )	p=b,b=c,c=p;
			if ( a>b )	p=a,a=b,b=p;
			tris[ntri].vertices[0]=points+a;
			tris[ntri].vertices[1]=points+b;
			tris[ntri].vertices[2]=points+c;
			ntri++;
		 }
	 }
  }

 qsort(tris,ntri,sizeof(triangle),compare_triangle_vertices);

 wtris=(triangle *)malloc(sizeof(triangle)*atri);
 memcpy(&wtris[0],&tris[0],sizeof(triangle));
 nwtri=1;

 for ( i=1 ; i<ntri ; i++ )
  {	if ( memcmp(tris[i-1].vertices,tris[i].vertices,3*sizeof(tpoint*)) )
	 {	memcpy(&wtris[nwtri],&tris[i],sizeof(triangle));
		nwtri++;
	 }
  }

 if ( rtris != NULL )	*rtris=wtris;
 if ( rntri != NULL )	*rntri=nwtri;

 free(tris);

 return(0);
}

/*****************************************************************************/
                                                 
