/*****************************************************************************/
/* maskdraw.c								     */
/*****************************************************************************/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <ctype.h>
#include <stdarg.h>

#include <fits/fits.h>

#include "io/tokenize.h"
#include "fitsmask.h"
#include "fi.h"

/*****************************************************************************/

#define		MASKDRAW_PIXEL		1
#define		MASKDRAW_BLOCK		2
#define		MASKDRAW_LINE		3
#define		MASKDRAW_CIRCLE		4

/*****************************************************************************/

static int maskdraw_draw_block(char **mask,int sx,int sy,
	int maskval,int x0,int y0,int bx,int by)
{
 int	k,l;

 for ( k=0 ; k<by ; k++ )
  {	if ( y0+k < 0 || y0+k >= sy )
		continue;
	for ( l=0 ; l<bx ; l++ )
	 {	if ( x0+l < 0 || x0+l >= sx )
			continue;
		mask[y0+k][x0+l] |= maskval;
	 }
  }
 return(0); 
}

static int maskdraw_draw_pixel(char **mask,int sx,int sy,
	int maskval,int x0,int y0)
{
 maskdraw_draw_block(mask,sx,sy,maskval,x0,y0,1,1);
 return(0);
}

static int maskdraw_draw_line(char **mask,int sx,int sy,
	int maskval,int x1,int y1,int x2,int y2,int width)
{
 int	len,lx,ly,i,x,y;

 lx=(x2>x1?x2-x1:x1-x2); 
 ly=(y2>y1?y2-y1:y1-y2); 
 len=(lx>ly?lx:ly);

 if ( width<0 )	width=0;

 if ( len<=0 )
  {	maskdraw_draw_block(mask,sx,sy,maskval,x1-width,y1-width,2*width+1,2*width+1);
	return(0);
  }

 for ( i=0 ; i<=len ; i++ )
  {	x=((2*x1+1)+2*(x2-x1)*i/len)/2;
	y=((2*y1+1)+2*(y2-y1)*i/len)/2;
	maskdraw_draw_block(mask,sx,sy,maskval,x-width,y-width,2*width+1,2*width+1);
  }

 return(0);
}

static int maskdraw_draw_circle(char **mask,int sx,int sy,
	int maskval,int x1,int y1,int radius)
{
 int	k,l,dx,dy2,r2;

 r2=radius*(radius+1);

 for ( k=y1-radius ; k<=y1+radius ; k++ )
  {	if ( ! ( 0<=k && k<sy ) )
		continue;
	dy2=(k-y1)*(k-y1);
	for ( l=x1-radius ; l<=x1+radius ; l++ )
	 { 	if ( ! ( 0<=l && l<sx ) )
			continue;
		dx=l-x1;
		if ( dy2+dx*dx<=r2 )
			mask[k][l] |= maskval;
	 }
  }
 return(0); 
 
}


/*****************************************************************************/

int maskdraw_parse_and_draw(char **mask,int sx,int sy,char *maskblock)
{
 char	*mtmp,*cmd[8],**cc;
 int	n,is_invalid,maskval,x0,y0,bx,by,br,w,type;

 if ( maskblock == NULL )
	return(1);
 else if ( mask == NULL || sx<=0 || sy<=0 )
	return(-1);

 mtmp=strdup(maskblock);
 n=tokenize_char(mtmp,cmd,':',7);
 is_invalid=0;

 x0=y0=0;
 bx=by=-1;
 maskval=0;

 if ( strcmp(cmd[0],"block")==0 )
	type=MASKDRAW_BLOCK,cc=cmd+1,n--;
 else if ( strcmp(cmd[0],"pixel")==0 )
	type=MASKDRAW_PIXEL,cc=cmd+1,n--;
 else if ( strcmp(cmd[0],"line")==0 )
	type=MASKDRAW_LINE,cc=cmd+1,n--;
 else if ( strcmp(cmd[0],"circle")==0 )
	type=MASKDRAW_CIRCLE,cc=cmd+1,n--;
 else 
	type=MASKDRAW_BLOCK,cc=cmd;

 if ( type==MASKDRAW_BLOCK )
  {	if ( n<2 || n>3 )
		is_invalid=1;
	else if ( (maskval=parse_mask_flags(cc[0]))<0 )
		is_invalid=2;
	else if ( sscanf(cc[1],"%d,%d",&x0,&y0)<2 )
	 	is_invalid=3;
	else if ( n>=3 && sscanf(cc[2],"%d,%d",&bx,&by)<2 )
		is_invalid=4;
	else
		is_invalid=0;

	if ( bx>=0 && by>=0 )
	 {	bx=bx-x0+1;
		by=by-y0+1;
	 }
	else
	 {	bx=1;
		by=1;
	 }

	if ( ! is_invalid )
		maskdraw_draw_block(mask,sx,sy,maskval,x0,y0,bx,by);
  }
 else if ( type==MASKDRAW_PIXEL )
  {	if ( n != 2  )
		is_invalid=1;
	else if ( (maskval=parse_mask_flags(cc[0]))<0 )
		is_invalid=2;
	else if ( sscanf(cc[1],"%d,%d",&x0,&y0)<2 )
	 	is_invalid=3;
	else
		is_invalid=0;

	if ( ! is_invalid )
		maskdraw_draw_pixel(mask,sx,sy,maskval,x0,y0);
  }
 else if ( type==MASKDRAW_LINE )
  {	w=0;

	if ( n<3 || n>4  )
		is_invalid=1;
	else if ( (maskval=parse_mask_flags(cc[0]))<0 )
		is_invalid=2;
	else if ( sscanf(cc[1],"%d,%d",&x0,&y0)<2 )
	 	is_invalid=3;
	else if ( sscanf(cc[2],"%d,%d",&bx,&by)<2 )
	 	is_invalid=3;
	else if ( n==4 && ( sscanf(cc[3],"%d",&w)<1 || w<0 ) )
	 	is_invalid=3;
	else
		is_invalid=0;

	if ( ! is_invalid )
		maskdraw_draw_line(mask,sx,sy,maskval,x0,y0,bx,by,w);
  }
 else if ( type==MASKDRAW_CIRCLE )
  {	w=0;

	if ( n != 3 )
		is_invalid=1;
	else if ( (maskval=parse_mask_flags(cc[0]))<0 )
		is_invalid=2;
	else if ( sscanf(cc[1],"%d,%d",&x0,&y0)<2 )
	 	is_invalid=3;
	else if ( sscanf(cc[2],"%d",&br)<1 )
	 	is_invalid=3;
	else
		is_invalid=0;

	if ( ! is_invalid )
		maskdraw_draw_circle(mask,sx,sy,maskval,x0,y0,br);
  }
 else
	is_invalid=5;

 free(mtmp);

 return(is_invalid);
}

/*****************************************************************************/
                 
