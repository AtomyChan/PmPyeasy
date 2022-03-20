/*****************************************************************************/
/* fitsmask.c 								     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Mask manipulation extensions to the 'fits.a' library.		     */
/* (c) 2004-05, Pal, A. (apal@szofi.elte.hu).				     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* See function prototypes and the usage of these and other related 	     */
/* functions in fitsmask.h and fits.h. The library is not standalone,	     */
/* it depends on fits.a.						     */
/*****************************************************************************/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include <fits/fits.h>

#include "fitsmask.h"

/*****************************************************************************/

static	char * maskhdr_default="MASKINFO";
#define		FITS_MASK_MAX	(0x7F)
#define		FITS_MASK_DEF	(0x01)

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/*
char **fits_mask_alloc(int sx,int sy)
{
 int	tsize,psize,bsize,cd,cb,snext,pnext;
 int	i,j,rank,arr[2],typesize;
 void	*ret,**pret;

 rank=2;
 typesize=sizeof(char);
 arr[0]=sx;
 arr[1]=sy;

 psize=0;bsize=1;
 for ( i=rank-1 ; i>=1 ; i-- )
  {	bsize=bsize*arr[i];
	psize+=bsize;
  }
 psize=psize*sizeof(void *);
 tsize=psize+typesize*bsize*arr[0];
 pret=(void **)malloc(tsize);
 if ( pret==NULL )	return(NULL);

 psize=0;bsize=1;
 for ( i=rank-1 ; i>=1 ; i-- )
  {	cd=arr[i];
	if ( i>1 )	snext=sizeof(void *);
	else		snext=typesize;
	cb=cd*bsize;
	pnext=psize+cb;
	for ( j=0 ; j<cb ; j++ )
		pret[psize+j]=(void *)( ((char *)(&pret[pnext]))+j*snext*arr[i-1] );
	psize=pnext;
	bsize=cb;
  }
 ret=(void *)pret;
 return((char **)(ret));
}
*/

char ** fits_mask_alloc(int sx,int sy)
{
 char	**ret,*r0;
 int	i;

 ret=(char **)malloc(sizeof(char *)*(size_t)sy+(size_t)sx*(size_t)sy);
 r0=(char *)(ret+(size_t)sy);
 for ( i=0 ; i<sy ; i++ )
  {	ret[i]=r0+(size_t)i*(size_t)sx;			}
 
 return(ret);
}

char ** fits_mask_create_empty(int sx,int sy)
{
 char	**mask;
 int	i;	

 if ( sx<=0 || sy<=0 )	return(NULL);
 mask=fits_mask_alloc(sx,sy);
 if ( mask==NULL )	return(NULL);
 for ( i=0 ; i<sy ; i++ )
  {	memset(mask[i],0,sx);		}
 return(mask);
}
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

static int fits_mask_mask_from_header_native(char **mask,fitsheaderset *header,int sx,int sy,int y0,int y1,char *maskhdr)
{
 char		buff[80],*wp,*sp;
 int		i,l,nhdr,n,x,y,lx,ly,xprev,yprev,use_diff,data;
 fitsheader	*fh;

 if ( maskhdr==NULL )	maskhdr=maskhdr_default;

 if ( sx<=0 || sy<=0 )	return(1);

 xprev=yprev=0;use_diff=0;

 nhdr=fits_headerset_get_count(header,maskhdr);
 data=FITS_MASK_DEF;
 for ( i=0 ; i<nhdr ; i++ )
  {	fh=fits_headerset_get_header(header,maskhdr,i);
	if ( fh->vtype != FITS_VSTR )	continue;
	strcpy(buff,fh->vstr);
	wp=buff;
	while ( *wp )
	 {	n=sscanf(wp,"%d,%d:%d,%d",&x,&y,&lx,&ly);
		if ( n==0 )		x=y=lx=ly=0;
		else if ( n==1 )
		 {	if ( x>0 )	use_diff=1;
			else if ( x<0 )	data=((-x)&FITS_MASK_MAX);
			else		use_diff=0;
			x=y=lx=ly=0;
		 }
		else if ( n==2 )	lx=ly=1;
		else if ( n==3 )
		 {	if ( lx>1 )		ly=1;
			else if ( lx<-1 )	ly=-lx,lx=1;
			else			lx=ly=1;
		 }

		if ( lx>0 && ly>0 )
		 {	if ( use_diff )		x+=xprev,y+=yprev;
			if ( x<0 )	lx+=x,x=0;
			if ( y<0 )	ly+=y,y=0;
			if ( x+lx>=sx )	lx=sx-x;
			if ( y+ly>=sy )	ly=sy-y;
			xprev=x,yprev=y;
			for ( ; ly>0 && lx>0 ; y++,ly-- )
			 {	if ( y<y0 || y>=y1 )		continue;
				for ( sp=mask[y-y0]+x,l=lx ; l>0 ; l--,sp++ )
					(*sp) |= data;
			 }
		 }
		while ( *wp && *wp != 32 )	wp++;
		while ( *wp==32 )		wp++;
	 };
  }
 return(0);
}

int fits_mask_mask_from_header(char **mask,fitsheaderset *header,int sx,int sy,char *maskhdr)
{
 int	r;
 if ( sx<=0 || sy<=0 )	return(1);
 r=fits_mask_mask_from_header_native(mask,header,sx,sy,0,sy,maskhdr);
 return(r);
}

char **fits_mask_read_from_header(fitsheaderset *header,int sx,int sy,char *maskhdr)
{
 char	**mask;
 int	r;

 if ( sx<=0 || sy<=0 )	return(NULL);

 mask=fits_mask_create_empty(sx,sy);
 if ( mask==NULL )	return(NULL);

 r=fits_mask_mask_from_header(mask,header,sx,sy,maskhdr);
 if ( r )
  {	fits_mask_free(mask);
	return(NULL);
  }
 else
	return(mask);
}

int fits_mask_mask_line(char *line,fitsheaderset *header,int sx,int sy,int y0,char *maskhdr)
{
 int	r;
 r=fits_mask_mask_from_header_native(&line,header,sx,sy,y0,y0+1,maskhdr);
 return(r);
}

int fits_mask_mask_more_line(char **lines,fitsheaderset *header,int sx,int sy,int y0,int ny,char *maskhdr)
{
 int	r;
 r=fits_mask_mask_from_header_native(lines,header,sx,sy,y0,y0+ny,maskhdr);
 return(r);
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

typedef struct
 {	int	x,y;
	int	sx,sy;
	int	value;
 } block;

#define	BLOCK_FRAGMENTS	128

static int fits_mask_scan_horizontal(char **mask,int sx,int sy,block **rblocks,int *rnblock)
{
 block	*blocks,*wb;
 int	i,j,oj,nj,nblock,nablock,pmv;
 char	*line;
 
 nablock=BLOCK_FRAGMENTS;
 blocks=(block *)malloc(sizeof(block)*nablock);nblock=0;

 for ( i=0 ; i<sy ; i++ )
  {	line=mask[i];
	for ( j=0 ; j<sx ; )
	 {	while ( j<sx && (! line[j]) )	j++;
		oj=j;
		pmv=(line[j]&FITS_MASK_MAX);
		while ( j<sx && (line[j]&FITS_MASK_MAX)==pmv )	j++;
		nj=j-oj;
		if ( nj>0 )
		 {	if ( nblock>=nablock )
			 {	nablock+=BLOCK_FRAGMENTS;
				blocks=(block *)realloc(blocks,sizeof(block)*nablock);
			 }
			wb=&blocks[nblock];
			wb->x=oj,wb->y=i;
			wb->sx=nj,wb->sy=1;
			wb->value=pmv;
			nblock++;
		 }
	 }
  }
 *rblocks=blocks,
 *rnblock=nblock;

 return(0);
}

/** old join_blocks algorithm... slow, but seems ok...			    **/
/*static int fits_mask_join_blocks(block *blocks,int nblock)
{
 int	i,j,fi,x,y,sx;
 for ( i=1 ; i<nblock ; i++ )
  {	sx=blocks[i].sx;
	if ( sx<=0 )	continue;
	x=blocks[i].x,y=blocks[i].y;
	for ( j=0,fi=-1 ; j<i && fi<0 ; j++ )
	 {	if ( blocks[j].x != x || blocks[j].sx != sx )	continue;
		if ( y == blocks[j].y+blocks[j].sy ) fi=j;
	 }
	if ( fi>=0 )
	 {	blocks[fi].sy+=blocks[i].sy;
		blocks[i].sx=-1;
	 }
  }
 return(0);	
}*/

/** new join_blocks algorithm... faster, but probably buggy...		    **/
static int fits_mask_join_blocks(block *blocks,int nblock)
{
 int	i,j,x,y,sx,sy,v;
 for ( i=0 ; i<nblock ; i++ )
  {	sx=blocks[i].sx;
	if ( sx<=0 )	continue;
	x=blocks[i].x,y=blocks[i].y;
	v=blocks[i].value;
	sy=blocks[i].sy;
	for ( j=i+1 ; j<nblock ; j++ )
	 {	if ( blocks[j].y < y+sy )	continue;
		else if ( blocks[j].y==y+sy )
		 {	if ( blocks[j].x==x && blocks[j].sx==sx && blocks[j].value==v )
			 {	blocks[i].sy+=blocks[j].sy;
				sy=blocks[i].sy;
				blocks[j].sx=-1;
			 }
			else	continue;
		 }
		else	break;
	 }
  }
 return(0);	
}

int fits_mask_export_as_header(fitsheaderset *header,int clean,char **mask,int sx,int sy,char *maskhdr)
{
 int	i,l,mxl,nblock,len,use_diff,x,y,xprev,yprev,vprev;
 block	*blocks,*wb;
 char	hdrstr[80],buff[80];

 if ( maskhdr==NULL )	maskhdr=maskhdr_default;
 if ( mask==NULL )	return(0);

 use_diff=1;

 if ( sx<=0 || sy<=0 )	return(1);

 if ( clean )
  {	fits_headerset_delete_all(header,maskhdr);		}
 
 fits_mask_scan_horizontal(mask,sx,sy,&blocks,&nblock);
 fits_mask_join_blocks(blocks,nblock); 
 sprintf(hdrstr,"%d",use_diff);len=strlen(hdrstr);
 mxl=67;xprev=yprev=0;vprev=-1;

/* for ( i=0,wb=blocks ; i<nblock ; i++,wb++ )
  {	fprintf(stderr,"[%d,%d] %dx%d (%d)\n",wb->x,wb->y,wb->sx,wb->sy,wb->value);	}
*/

 for ( i=0,wb=blocks ; i<nblock ; i++,wb++ )
  {	if ( wb->sx <= 0 )	continue;
	if ( ! use_diff )	x=wb->x,y=wb->y;
	else			x=wb->x-xprev,y=wb->y-yprev;

	if ( vprev != wb->value )
	 {	l=sprintf(buff,"-%d ",wb->value);
		vprev=wb->value;
	 }
	else
		l=0;

	if ( wb->sx==1 && wb->sy==1 )
		l+=sprintf(buff+l,"%d,%d",x,y);
	else if ( wb->sx>1 && wb->sy==1 )
		l+=sprintf(buff+l,"%d,%d:%d",x,y,wb->sx);
	else if ( wb->sx==1 && wb->sy>1 )
		l+=sprintf(buff+l,"%d,%d:%d",x,y,-wb->sy);
	else
		l+=sprintf(buff+l,"%d,%d:%d,%d",x,y,wb->sx,wb->sy);

	xprev=wb->x,yprev=wb->y;
	if ( len==0 )
	 {	strcpy(hdrstr,buff);
		len=l;
	 }
	else if ( len+l+1 <= mxl )
	 {	hdrstr[len]=32;strcpy(hdrstr+len+1,buff);
		len+=l+1;
	 }
	else
	 {	fits_headerset_set_string(header,maskhdr,FITS_SH_ADD,hdrstr,NULL);
		strcpy(hdrstr,buff);
		len=l;
	 }
  }
 if ( len>0 )
	fits_headerset_set_string(header,maskhdr,FITS_SH_ADD,hdrstr,NULL);

 free(blocks);

 return(0);
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

char ** fits_mask_expand_false(char **imask,int sx,int sy,int dis,int imv,int smv,int expand_border)
{
 int	i,j,k,l;
 char	**omask;

 omask=fits_mask_alloc(sx,sy);
 if ( omask==NULL )	return(NULL);
 
 for ( i=0 ; i<sy ; i++ )
  {	memcpy(omask[i],imask[i],sx);		}

 for ( i=0 ; i<sy ; i++ )
  {	for ( j=0 ; j<sx ; j++ )
	 {	if ( expand_border && ( i<dis || j<dis || i>=sy-dis || j>=sx-dis ) )
			omask[i][j] |= smv;
		else if ( imask[i][j] & imv )
		 {	for ( k=-dis ; k<=dis ; k++ )
			 {  if ( i+k<0 || i+k>=sy ) continue;
			    for ( l=-dis ; l<=dis ; l++ )
			     {	if ( j+l<0 || j+l>=sx ) continue;
				omask[i+k][j+l] |= smv;
			     }
			 }
		 }
	 }
  }
 return(omask);
}

int fits_mask_expand_logic(char **imask,int sx,int sy,int dis,int imv,int smv,int expand_border)
{
 int	i,j,k,l;
 char	**omask;

 omask=fits_mask_alloc(sx,sy);
 if ( omask==NULL )	return(-1);
 
 for ( i=0 ; i<sy ; i++ )
  {	memcpy(omask[i],imask[i],sx);		}

 for ( i=0 ; i<sy ; i++ )
  {	for ( j=0 ; j<sx ; j++ )
	 {	if ( expand_border && ( i<dis || j<dis || i>=sy-dis || j>=sx-dis ) )
			omask[i][j] |= expand_border;
		else if ( imask[i][j] & imv )
		 {	int	ival;
			ival = (imask[i][j] & imv) | smv;
			for ( k=-dis ; k<=dis ; k++ )
			 {  if ( i+k<0 || i+k>=sy ) continue;
			    for ( l=-dis ; l<=dis ; l++ )
			     {	if ( j+l<0 || j+l>=sx ) continue;
				omask[i+k][j+l] |= ival;
			     }
			 }
		 }
	 }
  }

 for ( i=0 ; i<sy ; i++ )
  {	memcpy(imask[i],omask[i],sx);		}

 fits_mask_free(omask);

 return(0);
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

char **fits_mask_duplicate(char **mask,int sx,int sy)
{
 char 	**nm,**wm;
 nm=fits_mask_alloc(sx,sy);
 if ( nm==NULL )	return(NULL);
 for ( wm=nm ; sy>0 ; sy--,wm++,mask++ )
  {	memcpy(*wm,*mask,sx);			}
 return(nm);
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

int fits_mask_and(char **mask,int sx,int sy,char **other)
{
 char	*p,*q;
 int	s;
 if ( other==NULL )	return(0);
 for ( ; sy>0 ; sy--,mask++,other++ )
  {	for ( p=*mask,q=*other,s=sx ; s>0 ; s--,p++,q++ )
	 {	(*p) |= (*q);				}
  }
 return(0);
}

int fits_mask_or(char **mask,int sx,int sy,char **other)
{
 char	*p,*q;
 int	s;
 if ( other==NULL )	return(0);
 for ( ; sy>0 ; sy--,mask++,other++ )
  {	for ( p=*mask,q=*other,s=sx ; s>0 ; s--,p++,q++ )
	 {	(*p) = ~((~(*p))|(~(*q)));			}
  }
 return(0);
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

char ** fits_mask_create_floyd(int sx,int sy,int a,int b,int svm)
{
 char 	**mask;
 int	i,j,r,sgn,n,w,c,*im[2],*impnt,max,*wp,*l0,*l1,offset;

 if ( b<=0 )		return(NULL);
 if ( a<0 || a>b )	return(NULL);

 if ( sx<=0 || sy<=0 )	return(NULL); 

 mask=fits_mask_create_empty(sx,sy);
 if ( mask==NULL )	return(NULL);
 impnt=(int *)malloc(sizeof(int)*sx*2);
 im[0]=impnt,im[1]=impnt+sx;

 if ( impnt==NULL )
  {	fits_mask_free(mask);return(NULL);	}

 max=32768;l0=im[0],l1=im[1];
 for ( j=0,r=a*max ; j<sx ; j++ )
  {	l0[j]=r/b;
	r=r%b+a*max;
  }
 offset=sx;
 if ( a>0 )	offset+=2*b/a;
 for ( i=-offset,sgn=1 ; i<sy ; i++,sgn=-sgn )
  {	for ( j=0 ; j<sx ; j++ )
	 {	l1[j]=r/b;
		r=r%b+a*max;
	 }
	if ( sgn>0 )	
	 {	for ( j=0 ; j<sx ; j++ )
		 {	n=5;
			if ( j<sx-1 )	n+=8;
			if ( j>0 )	n+=3;
			if ( 2*l0[j]<max )	c=0;
			else			c=max;
			w=l0[j]-c;
			l1[j]+=5*w/n;
			if ( j<sx-1 )	l0[j+1]+=7*w/n,l1[j+1]+=w/n;
			if ( j>0 )	l1[j-1]+=3*w/n;
		 }
	 }
	else
	 {	for ( j=sx-1 ; j>=0 ; j-- )
		 {	n=5;
			if ( j>0 )	n+=8;
			if ( j<sx-1 )	n+=3;
			if ( 2*l0[j]<max )	c=0;
			else			c=max;
			w=l0[j]-c;
			l1[j]+=5*w/n;
			if ( j>0 )	l0[j-1]+=7*w/n,l1[j-1]+=w/n;
			if ( j<sx-1 )	l1[j+1]+=3*w/n;
		 }
	 }
	if ( i>=0 )
	 {	for ( j=0 ; j<sx ; j++ )
		 {	if ( 2*l0[j]>=max )	mask[i][j]=0;
			else			mask[i][j]=svm;
		 }
	 }
	wp=l0,l0=l1,l1=wp;
  }
 free(impnt);

 return(mask);
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

int fits_mask_mark_nans(fitsimage *img,char **mask,int setmask)
{
 int	i,j,sx,sy;
 if ( mask==NULL )			return(-1);
 if ( img==NULL || img->data==NULL )	return(-1);
 sx=img->sx,sy=img->sy;
 if ( sx<=0 || sy<=0 )			return(-1);

 for ( i=0 ; i<sy ; i++ )
  {	for ( j=0 ; j<sx ; j++ )
	 {	if ( ! isfinite(img->data[i][j]) )	
		 {	mask[i][j] |= setmask;
			img->data[i][j]=0.0;
		 }
	 }
  }

 return(0);
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

static int fits_mask_expand(char ***rmask,int sx,int sy,int nx,int ny,int fill)
{
 char	**mask;
 size_t	ns;
 int	i;

 if ( nx<sx || ny<sy )
	return(-1);
 else if ( nx==sx && ny==sy )
	return(0);

 mask=*rmask;

 ns=sizeof(char *)*(size_t)ny+(size_t)nx*(size_t)ny;
 mask=realloc(mask,ns);

 for ( i=sy-1 ; i>=0 ; i-- )
  {	memmove(((char *)(mask+ny))+i*nx,((char *)(mask+sy))+i*sx,sx);
	if ( nx>sx )	memset(((char *)(mask+ny))+i*nx+sx,fill,nx-sx);
  }
 for ( i=sy ; i<ny ; i++ )
	memset(((char *)(mask+ny))+i*nx,fill,nx);

 for ( i=0 ; i<ny ; i++ )
	mask[i]=(char *)(mask+ny)+i*nx;

 *rmask=mask;

 return(0); 
}

static int fits_mask_shrink(char ***rmask,int sx,int sy,int nx,int ny)
{
 char	**mask;
 size_t	ns;
 int	i;

 if ( nx>sx || ny>sy )
	return(-1);

 else if ( nx==sx && ny==sy )
	return(0);

 mask=*rmask;

 for ( i=0 ; i<ny ; i++ )
  {	memmove(((char *)(mask+ny))+i*nx,((char *)(mask+sy))+i*sx,nx);
	mask[i]=(char *)(mask+ny)+i*nx;
  }

 ns=sizeof(char *)*(size_t)ny+(size_t)nx*(size_t)ny;
 mask=realloc(mask,ns);

 *rmask=mask;
 
 return(0);
}

int fits_mask_copy_line(char **mask,int tx,int ty,int x0,int y0,char *line,int outer)
{
 if ( y0<0 || y0>=ty )
  {	memset(line,outer,tx);
	return(0);
  }
 else
  {	int	l,c;
	l=tx;
	while ( l>0 )
	 {	if ( x0<0 )
		 {	c=-x0;
			if ( c>l )	c=l;
			memset(line,outer,c);
			line+=c;
			l-=c;
			x0+=c;
		 }
		else if ( x0>=tx )
		 {	c=l;
			memset(line,outer,c);
			line+=c;
			l-=c;
			x0+=c;
		 }
		else
		 {	c=tx-x0;
			if ( c>l )	c=l;
			memcpy(line,mask[y0]+x0,c);
			line+=c;
			l-=c;
			x0+=c;
		 }
	 }
	return(0);
  }
}

int fits_mask_trim(char ***rmask,int sx,int sy,int x0,int y0,int nx,int ny,int outer)
{
 char	**mask,*line;
 int	tx,ty,i;

 mask=*rmask;

 tx=(nx>sx?nx:sx);
 ty=(ny>sy?ny:sy);

 fits_mask_expand(&mask,sx,sy,tx,ty,outer);

 line=(char *)malloc(tx);

 if ( x0 != 0 || y0 != 0 )
  {	if ( y0>0 )
	 {	for ( i=0 ; i<ty ; i++ )
		 {	fits_mask_copy_line(mask,tx,ty,x0,i+y0,line,outer);
			memcpy(mask[i],line,tx);
		 }
	 }
	else
	 {	for ( i=ty-1 ; i>=0 ; i-- )
		 {	fits_mask_copy_line(mask,tx,ty,x0,i+y0,line,outer);
			memcpy(mask[i],line,tx);
		 }
	 }
  }

 free(line);

 fits_mask_shrink(&mask,tx,ty,nx,ny); 

 *rmask=mask;
 
 return(0);
}

/*****************************************************************************/

int fits_mask_is_clean(char **mask,int sx,int sy)
{ 
 char	*line;
 int	i;

 line=(char *)malloc(sx);
 memset(line,0,sx);

 for ( i=0 ; i<sy ; i++ )
  {	if ( memcmp(mask[i],line,sx) )
	 {	free(line);
		return(0);
	 }
  }
 free(line);
 return(1);
}

/*****************************************************************************/

int fits_mask_free(char **mask)
{
 if ( mask != NULL )	free(mask);
 return(0);
}

/*****************************************************************************/
                                                             
