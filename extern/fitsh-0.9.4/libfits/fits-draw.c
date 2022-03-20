/*****************************************************************************/
/* fits-draw.c 								     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Another simple standalone library for manipulating FITS files:	     */
/* Drawing pixels, circles and lines to 2-dimensional images.		     */ 
/* (c) 2004-06, Pal, A. (apal@szofi.elte.hu). 				     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* See function prototypes and the usage of the functions in fits.h          */
/*****************************************************************************/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>

#define	_FITS_SOURCE

#include <fits/fits.h>
#include "fits-common.h"

/*****************************************************************************/

int fits_image_draw_pixel(fitsimage *img,int x,int y,double c)
{
 if ( img==NULL || img->vdata==NULL )		return(1);
 if ( img->data==NULL || img->dim != 2 )	return(1);
 
 if ( x<0 || y<0 || x>=img->sx || y>=img->sy )	return(0);
 img->data[y][x]=c;
 return(0);
}

int fits_image_draw_pixel_nochk(fitsimage *img,int x,int y,double c)
{
 if ( x<0 || y<0 || x>=img->sx || y>=img->sy )	return(0);
 img->data[y][x]=c;
 return(0);
}

int fits_image_draw_line(fitsimage *img,int x,int y,int dx,int dy,double col,int style)
{
 int sx,sy,s,n;

 if ( img==NULL || img->vdata==NULL )		return(1);
 if ( img->data==NULL || img->dim != 2 )	return(1);

 style = style & 0xffff;
 style = style | (style<<16);

 if ( style & 1 )	fits_image_draw_pixel_nochk(img,x,y,col);

 dx-=x,dy-=y;
 if ( dx==0 && dy==0 )	return(0);
 
 if ( dx<0 )	dx=-dx,sx=-1;	else	sx=+1;
 if ( dy<0 )	dy=-dy,sy=-1;	else	sy=+1;

 if ( dx>=dy )
  {	for ( s=dx/2,n=dx ; n ; n-- )
	 {	style=(style<<1)|((style>>31)&1);
		if ( style & 1 )	fits_image_draw_pixel_nochk(img,x,y,col);
		s+=dy;
		if ( s>=dx )	y+=sy,s-=dx;
		x+=sx;
	 };
  }
 else
  {	for ( s=dy/2,n=dy ; n ; n-- )
	 {	style=(style<<1)|((style>>31)&1);
		if ( style & 1 )	fits_image_draw_pixel_nochk(img,x,y,col);
		s+=dx;
		if ( s>=dy )	x+=sx,s-=dy;
		y+=sy;
	 };
  }	
 return(0);
}

int fits_image_draw_circle(fitsimage *img,int x,int y,int r,double col)
{
 int	a,b,c,d,dx,dy;

 if ( img==NULL || img->vdata==NULL )		return(1);
 if ( img->data==NULL || img->dim != 2 )	return(1);

 dx=r,dy=0;
 a=2*r-1,b=1,c=3,d=2*r-3;

 while ( dx>=dy )
  {	fits_image_draw_pixel_nochk(img,x+dx,y+dy,col);
	fits_image_draw_pixel_nochk(img,x+dx,y-dy,col);
 	fits_image_draw_pixel_nochk(img,x-dx,y+dy,col);
	fits_image_draw_pixel_nochk(img,x-dx,y-dy,col);
	if ( dx != dy )
	 {	fits_image_draw_pixel_nochk(img,x+dy,y+dx,col);
		fits_image_draw_pixel_nochk(img,x+dy,y-dx,col);
 		fits_image_draw_pixel_nochk(img,x-dy,y+dx,col);
		fits_image_draw_pixel_nochk(img,x-dy,y-dx,col);
	 }
	if ( a>=b )	b+=c,c+=2,dy++;
	if ( a<=b )	a+=d,d-=2,dx--;
  };

 return(0);
}

/*****************************************************************************/

int fits_draw_pixel(fits *img,int x,int y,double c)
{
 int	ret;
 ret=fits_image_draw_pixel(&img->i,x,y,c);
 return(ret);
}

int fits_draw_line(fits *img,int x,int y,int dx,int dy,double col,int style)
{
 int	ret;
 ret=fits_image_draw_line(&img->i,x,y,dx,dy,col,style);
 return(ret);
}

int fits_draw_circle(fits *img,int x,int y,int r,double col)
{
 int	ret;
 ret=fits_image_draw_circle(&img->i,x,y,r,col);
 return(ret);
}

/*****************************************************************************/
                                                        
                   
