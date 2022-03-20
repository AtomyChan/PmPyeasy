/*****************************************************************************/
/* fits-header.c 							     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Another simple standalone library for manipulating FITS files:	     */
/* FITS header manipulation.					             */
/* (c) 2004-06, Pal, A. (apal@szofi.elte.hu). 				     */
/* See reference(s) at the end of this source code.			     */
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

int fits_headerset_reset(fitsheaderset *header)
{
 header->hdrs=NULL;
 header->nhdr=0;
 header->ahdr=0;
 return(0);
}
int fits_headerset_duplicate(fitsheaderset *ret,fitsheaderset *act)
{
 ret->nhdr=act->nhdr;
 ret->ahdr=act->ahdr;
 ret->hdrs=(fitsheader *)malloc(sizeof(fitsheader)*(act->ahdr));
 memcpy(ret->hdrs,act->hdrs,sizeof(fitsheader)*(act->nhdr));
 return(0);
}
int fits_headerset_free(fitsheaderset *header)
{
 if ( header->hdrs != NULL )
  {	free(header->hdrs);		}
 fits_headerset_reset(header);
 return(0);
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

static fitsheader *fits_headerset_add(fitsheaderset *header)
{ 
 fitsheader	*fhd;

 if ( header->hdrs==NULL || header->nhdr==0 || header->ahdr==0 )
  {	header->hdrs=(fitsheader *)malloc(sizeof(fitsheader)*HDRBLOCKS);
	header->nhdr=0;
	header->ahdr=HDRBLOCKS;
  }
 else if ( header->nhdr+1 > header->ahdr )
  {	header->hdrs=(fitsheader *)realloc(header->hdrs,sizeof(fitsheader)*(header->ahdr+HDRBLOCKS));
	header->ahdr+=HDRBLOCKS;
  }
 fhd=&header->hdrs[header->nhdr];
 memset(fhd,0,sizeof(fitsheader));
 header->nhdr++;

 return(fhd);
}

static int fits_headerset_import_values(char *buff,int i0,char *comment,int *rtype,int *rvint,char *vstr,double *rvdouble)
{
 int	i,j,k,l,type,vint;
 double	vdouble;

 type=-1;
 vint=0;
 vdouble=0.0;

 j=0;
 for ( i=i0 ; i<FITS_TAPE_CARDIMAGESIZE ; i++ )
  {	if ( buff[i]=='\'' )		j=!j;
	if ( buff[i]=='/' && !j )	break;
  }
 if ( i<FITS_TAPE_CARDIMAGESIZE )
  {	buff[i]=0;i++;
	while ( buff[i]==32 )	i++;
	strcpy(comment,buff+i);
	while ( strlen(comment)>0 && comment[strlen(comment)-1]==32 )
		comment[strlen(comment)-1]=0;
  }

 j=0;
 for ( i=i0 ; i<strlen(buff) ; i++ )
  {	if ( buff[i]=='-' || ('0'<=buff[i] && buff[i]<='9') || buff[i]=='.' || buff[i]=='+' )
 	 {	j=1;break;		}
	else if ( buff[i] != 32 )
 	 {	j=2;break;		}
  }

 /** Numerical card image value:		**/
 if ( j==1 )
  {	for ( k=i,l=0 ; buff[k] ; k++ )
	 {	if ( buff[k]=='.' )	l=1;	}
	sscanf(buff+i,"%lg",&vdouble);
	vint=(int)vdouble;
	if ( (double)vint==vdouble && ! l )
	 {	type=FITS_VINT;		}
	else
		type=FITS_VDOUBLE;
  }

 /** Otherwise, some kind of string or boolean:	**/
 else if ( j==2 )
  {	if ( buff[i]=='\'' )
	 {	type=FITS_VSTR;
		for ( j=i+1,k=0 ; buff[j] ; j++ )
		 {	if ( buff[j]=='\'' && buff[j+1]=='\'' )
			vstr[k]='\'',k++,j++;
			else if ( buff[j]=='\'' )
				break;
			else
				vstr[k]=buff[j],k++;
		 }
		vstr[k]=0;
	 }
	else
	 {	for ( j=i,k=0 ; buff[j] ; j++ )
			vstr[k]=buff[j],k++;
		vstr[k]=0;
		while ( k>0 && vstr[k-1]==32 )	k--,vstr[k]=0;
		if ( strcmp(vstr,"T")==0 )
			type=FITS_VBOOLEAN,vint=1;
		else if ( strcmp(vstr,"F")==0 )
			type=FITS_VBOOLEAN,vint=0;
		else
			type=FITS_VSTR;
	 }
   }

 *rtype = type;
 *rvint = vint;
 *rvdouble = vdouble;

 return(0);
}

int fits_headerset_read_cb(int (*cb_read)(void *,void *,int),void *param,fitsheaderset *header)
{
 char		buff[2*FITS_TAPE_CARDIMAGESIZE];
 char		name[32],comment[2*FITS_TAPE_CARDIMAGESIZE];
 int		i,k,lin,type;

 int		vint;
 double		vdouble;
 char		vstr[128];

 fitsheader	*fhd;

 lin=0;
 while ( 1 )
  {	/** Read the next 80 bytes from the stream. This block is	**/
	/** called "card image", see [1].				**/

	i=cb_read(param,buff,FITS_TAPE_CARDIMAGESIZE);
	if ( i <=0 )	break;

	lin++;
	for ( i=0,k=0 ; i<FITS_TAPE_CARDIMAGESIZE ; i++ )
	 {	if ( buff[i]==32 )	k++;		}

	/** Uncomment this if the end of the header can only be denoted	**/
	/** by the keyword 'END'. Otherwise end-of-header will occur	**/
	/** at the first totally empty header field. Note that according**/
	/** to FITS format, the keyword 'END' is mandatory (see [1]) !  **/
	if ( k==FITS_TAPE_CARDIMAGESIZE )
	 {	fhd=fits_headerset_add(header);
		fhd->vtype=FITS_EMPTY;
		continue;
	 }

	buff[FITS_TAPE_CARDIMAGESIZE]=0;

	/** The card image is absolutely empty or contains the reserved	**/
	/** keyword "END". Both cases indicate the end of the header.	**/
	if ( k==FITS_TAPE_CARDIMAGESIZE || memcmp(buff,"END     ",8)==0 )
	 {	while ( lin % (FITS_TAPE_CARDIMAGECOUNT) )
		 {	cb_read(param,buff,FITS_TAPE_CARDIMAGESIZE);lin++;	};
		break;
	 }

	/** Otherwise, its a normal card image with some context.	**/
	else
	 {	char	*p;

		/** First, the keyword should be read...		**/
		for ( i=0 ; i<8 && buff[i] != 32 ; i++ )
		 {	name[i]=buff[i];			}

		/** Reset the target storage variables...		**/ 
		name[i]=0;type=0;vint=0;vdouble=0.0;
		buff[FITS_TAPE_CARDIMAGESIZE]=0;comment[0]=0;
		
		/** Standard commentatory keywords: COMMENT and HISTORY **/
		/** Any ASCII text after the 8th col. is acceptable. 	**/
		if ( strcmp(name,"COMMENT")==0 || strcmp(name,"HISTORY")==0 )
		 {	strcpy(vstr,buff+8);
			type=FITS_VCOMMENT;
		 }

		else if ( strcmp(name,"HIERARCH")==0 && buff[8]==' ' && (p=strchr(buff,'=')) != NULL && *(p-1)==' ' && *(p+1)==' ' )
		 {	int	len,offset;
			len=(p-buff)-1-9;
			memcpy(name,buff+9,len);
			name[len]=0;
			offset=(p-buff)+2;
			fits_headerset_import_values(buff,offset,comment,&type,&vint,vstr,&vdouble);
		 }

		/** The card image has an appropriate value field,	**/
		/** indicated by the substring "= " at the position 8.	**/
		else if ( memcmp(buff+8,"= ",2)==0 )	
		 {	fits_headerset_import_values(buff,10,comment,&type,&vint,vstr,&vdouble);
		 }

		/** Any other case: the card image is a commentatory 	**/
		/** field without a value indicator.	     		**/
		else
		 {	strcpy(vstr,buff+8);
			type=FITS_VCOMMENT;
		 }

		if ( type<0 )
		 {	strcpy(vstr,buff+8);
			type=FITS_VCOMMENT;
		 }

		fhd=fits_headerset_add(header);
		
		fhd->vtype=type;

		if ( type==FITS_VINT || type==FITS_VBOOLEAN )
			fhd->vint=vint;
		else if ( type==FITS_VDOUBLE )
			fhd->vdouble=vdouble;
		else if ( type==FITS_VSTR || type==FITS_VCOMMENT )
			strcpy(fhd->vstr,vstr);

		strcpy(fhd->name,name);
		strcpy(fhd->comment,comment);
	 }
  };

 if ( header->hdrs != NULL )
  {	int	n;
	n=header->nhdr;
	while ( 0<n && header->hdrs[n-1].vtype==FITS_EMPTY )	n--;
	header->nhdr=n;
  }

 return(0);
}
int fits_headerset_read(FILE *fr,fitsheaderset *header)
{
 return(fits_headerset_read_cb(fits_cb_read,(void *)fr,header));
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

int fits_headerset_write_cb(int (*cb_write)(void *,void *,int),void *param,fitsheaderset *header)
{
 char	buff[256],val[256];
 int	i,j,k,align,ln,is_hierarch;

 ln=0;
 for ( i=0 ; i<header->nhdr ; i++ )
  {	
	align=0;
	switch ( header->hdrs[i].vtype )
	 {   case FITS_VINT:
		sprintf(val,"%d",header->hdrs[i].vint);
		align=2;
		break;
	     case FITS_VDOUBLE:
		sprintf(val,"%.15g", header->hdrs[i].vdouble);
		j=(int)header->hdrs[i].vdouble;
		if ( (double)j==header->hdrs[i].vdouble )
		 {	strcat(val,".");		}
		align=2;
		break;
	     case FITS_VSTR:
		if ( strcmp(header->hdrs[i].name,headers[HDR_XTENSION])==0 )
		 {	sprintf(val,"'%-8s'",header->hdrs[i].vstr);	}
		else
		 {	sprintf(val,"'%s'",header->hdrs[i].vstr);	}
		align=1;
		break;
	     case FITS_VBOOLEAN:
		if ( header->hdrs[i].vint )	strcpy(val,"T");
		else				strcpy(val,"F");
		align=2;
		break;
	     case FITS_VCOMMENT:
		strcpy(val,header->hdrs[i].vstr);
		align=3;
		break;
	     case FITS_EMPTY:
		align=-1;
		break;

	 }

	if ( ! align )	continue;

	if ( 8<strlen(header->hdrs[i].name) )	
		is_hierarch=1;
	else
		is_hierarch=0;

	if ( ! is_hierarch )
	 {
		strncpy(buff,header->hdrs[i].name,8);buff[8]=0;
		while ( strlen(buff)<8 )	strcat(buff," ");

		if ( align<0 )
		 {	memset(buff,32,FITS_TAPE_CARDIMAGESIZE);
			buff[FITS_TAPE_CARDIMAGESIZE]=0;
		 }

		else if ( align==1 )
		 {	strcat(buff,"= ");
			strcat(buff,val);
		 }
		else if ( align==2 )
		 {	strcat(buff,"= ");
			while ( strlen(buff)+strlen(val)<30 )
			 {	strcat(buff," ");		}
			strcat(buff,val);
		 }
		else if ( align==3 )
		 {	strcat(buff,val);			}

		if ( header->hdrs[i].comment[0] && ( align==1 || align==2 ) )
		 {	j=strlen(buff);
			k=strlen(header->hdrs[i].comment);
			while ( j<30 && j+k+3<FITS_TAPE_CARDIMAGESIZE )
			 {	buff[j]=32,buff[j+1]=0;
				j++;
			 };
			strcat(buff," / ");
			strcat(buff,header->hdrs[i].comment);
		 }

		buff[FITS_TAPE_CARDIMAGESIZE]=0;
		while ( strlen(buff)<FITS_TAPE_CARDIMAGESIZE )	strcat(buff," ");
	 }

	else
	 {	int	n;

		if ( header->hdrs[i].comment[0] )
		 {	n=snprintf(buff,FITS_TAPE_CARDIMAGESIZE,
				"HIERARCH %s = %s / %s",header->hdrs[i].name,val,header->hdrs[i].comment);
		 }
		else
		 {	n=snprintf(buff,FITS_TAPE_CARDIMAGESIZE+1,
				"HIERARCH %s = %s",header->hdrs[i].name,val);
		 }
		if ( n<FITS_TAPE_CARDIMAGESIZE )
			memset(buff+n,32,FITS_TAPE_CARDIMAGESIZE-n);

		buff[FITS_TAPE_CARDIMAGESIZE]=0;
	 }

	cb_write(param,buff,FITS_TAPE_CARDIMAGESIZE);
	ln++;
  }
 memset(buff,32,FITS_TAPE_CARDIMAGESIZE);
 memcpy(buff,headers[HDR_END],strlen(headers[HDR_END]));
 cb_write(param,buff,FITS_TAPE_CARDIMAGESIZE);
 ln++;

 memset(buff,32,FITS_TAPE_CARDIMAGESIZE);
 while ( (ln)%FITS_TAPE_CARDIMAGECOUNT )
  {	cb_write(param,buff,FITS_TAPE_CARDIMAGESIZE);
	ln++;
  };

 return(ln*FITS_TAPE_CARDIMAGESIZE);	
}
int fits_headerset_write(FILE *fw,fitsheaderset *header)
{
 return(fits_headerset_write_cb(fits_cb_write,(void *)fw,header));
}

/*****************************************************************************/

int fits_headerset_get_count(fitsheaderset *header,char *hdr)
{
 int	i,n; 
 if ( header==NULL || header->hdrs==NULL )	return(0);
 for ( i=0,n=0 ; i<header->nhdr ; i++ )
  {	if ( strcmp(header->hdrs[i].name,hdr)==0 )	n++;		}
 return(n);
}

int fits_headerset_get_id(fitsheaderset *header,char *hdr,int cnt)
{
 int	i,n;
 if ( header==NULL || header->hdrs==NULL )	return(-1);
 for ( i=0,n=0 ; i<header->nhdr ; i++ )
  {	if ( strcmp(header->hdrs[i].name,hdr)==0 )
	 {	if ( n==cnt )	return(i);
		else		n++;
	 }
  }
 return(-1);
}

fitsheader *fits_headerset_get_header(fitsheaderset *header,char *hdr,int cnt)
{
 int	i;
 i=fits_headerset_get_id(header,hdr,cnt);
 if ( i<0 )	return(NULL);
 else		return(&header->hdrs[i]);
}

fitsheader *fits_headerset_get_uniq_header(fitsheaderset *header,char *hdr)
{
 int	i,cnt;
 cnt=fits_headerset_get_count(header,hdr);
 if ( cnt != 1 )	return(NULL);
 i=fits_headerset_get_id(header,hdr,0);
 if ( i<0 )		return(NULL);
 else			return(&header->hdrs[i]);
}

int fits_headerset_get_as_double(fitsheaderset *header,char *hdr,double *ret,int is_ambigous_allowed)
{
 int		n;
 fitsheader	*h;
 if ( ! is_ambigous_allowed )
  {	n=fits_headerset_get_count(header,hdr);
	if ( n != 1 )	return(1);
  }
 h=fits_headerset_get_header(header,hdr,0);
 if ( h==NULL )	return(1);
 if ( h->vtype == FITS_VINT )
  {	*ret=(double)h->vint;
	return(0);
  }
 else if ( h->vtype == FITS_VDOUBLE )
  {	*ret=h->vdouble;
	return(0);
  }
 else
  {	return(1);	}
}


fitsheader *fits_headerset_append(fitsheaderset *img)
{
 fitsheader	*hdr;
 if ( img->hdrs==NULL || img->nhdr==0 || img->ahdr==0 )
  {	img->hdrs=(fitsheader *)malloc(sizeof(fitsheader)*HDRBLOCKS);
	img->nhdr=0;
	img->ahdr=HDRBLOCKS;
  }
 else if ( img->nhdr+1 > img->ahdr )
  {	img->hdrs=(fitsheader *)realloc(img->hdrs,sizeof(fitsheader)*(img->ahdr+HDRBLOCKS));
	img->ahdr+=HDRBLOCKS;
  }
 hdr=&img->hdrs[img->nhdr];
 img->nhdr++;
 return(hdr);
}

fitsheader * fits_headerset_insert_to(fitsheaderset *img,int n)
{
 if ( img->hdrs==NULL || img->nhdr==0 || img->ahdr==0 )
  {	img->hdrs=(fitsheader *)malloc(sizeof(fitsheader)*HDRBLOCKS);
	img->nhdr=1;
	img->ahdr=HDRBLOCKS;
	n=0;
  }
 else 
  {	if ( img->nhdr+1 > img->ahdr )
	 {	img->hdrs=(fitsheader *)realloc(img->hdrs,sizeof(fitsheader)*(img->ahdr+HDRBLOCKS));
		img->ahdr+=HDRBLOCKS;
	 }
	if ( n>img->nhdr )	n=img->nhdr;
	if ( n<img->nhdr )
	 {	memmove(img->hdrs+n+1,img->hdrs+n,sizeof(fitsheader)*(img->nhdr-n));	}
	img->nhdr++;
  }
 return(&img->hdrs[n]);
}

fitsheader * fits_headerset_insert(fitsheaderset *img)
{
 if ( img->hdrs==NULL || img->nhdr==0 || img->ahdr==0 )
  {	img->hdrs=(fitsheader *)malloc(sizeof(fitsheader)*HDRBLOCKS);
	img->nhdr=1;
	img->ahdr=HDRBLOCKS;
  }
 else 
  {	if ( img->nhdr+1 > img->ahdr )
	 {	img->hdrs=(fitsheader *)realloc(img->hdrs,sizeof(fitsheader)*(img->ahdr+HDRBLOCKS));
		img->ahdr+=HDRBLOCKS;
	 }
	memmove(img->hdrs+1,img->hdrs,sizeof(fitsheader)*(img->nhdr));
	img->nhdr++;
  }
 return(&img->hdrs[0]);
}

fitsheader *fits_headerset_set_any(fitsheaderset *img,char *hdr,int rule,char *comment)
{
 fitsheader	*hd;
 int		n;
 n=fits_headerset_get_count(img,hdr);
 if ( rule==FITS_SH_FORCEFIRST || rule==FITS_SH_INSERT )
  {	if ( rule==FITS_SH_FORCEFIRST )
		fits_headerset_delete_all(img,hdr);
	hd=fits_headerset_insert(img);
	strncpy(hd->name,hdr,15);
	hd->name[15]=0;
	hd->comment[0]=0;
	hd->vtype=0;
  }	
 else if ( rule==FITS_SH_BEGIN )
  {	for ( n=0 ; n<img->nhdr && img->hdrs != NULL ; n++ )
	 {	if ( strcmp(img->hdrs[n].name,headers[HDR_SIMPLE])==0 )
			continue;
		else if ( strcmp(img->hdrs[n].name,headers[HDR_XTENSION])==0 )
			continue;
		else if ( memcmp(img->hdrs[n].name,headers[HDR_NAXIS],5)==0 )
			continue;
		else if ( strcmp(img->hdrs[n].name,headers[HDR_BITPIX])==0 )
			continue;
		else if ( strcmp(img->hdrs[n].name,headers[HDR_EXTEND])==0 )
			continue;
		else	
			break;
	 }
	hd=fits_headerset_insert_to(img,n);
	strncpy(hd->name,hdr,15);
	hd->name[15]=0;
	hd->comment[0]=0;
	hd->vtype=0;
  }
 else if ( n==0 || rule==FITS_SH_ADD )
  {	hd=fits_headerset_append(img);
	strncpy(hd->name,hdr,15);
	hd->name[15]=0;
	hd->comment[0]=0;
	hd->vtype=0;
  }
 else
  {	if ( rule==FITS_SH_FIRST )	hd=fits_headerset_get_header(img,hdr,0);
	else if ( rule==FITS_SH_LAST )	hd=fits_headerset_get_header(img,hdr,n-1);
	else				hd=NULL;
  }

 if ( comment != NULL && hd != NULL )
  {	strncpy(hd->comment,comment,79);
	hd->comment[79]=0;
  }
 return(hd);
}

int fits_headerset_set_integer(fitsheaderset *img,char *hdr,int rule,int val,char *comment)
{
 fitsheader	*hd;
 hd=fits_headerset_set_any(img,hdr,rule,comment);
 hd->vtype=FITS_VINT;
 hd->vint=val;	
 return(0);
}
int fits_headerset_set_double(fitsheaderset *img,char *hdr,int rule,double val,char *comment)
{
 fitsheader	*hd;
 hd=fits_headerset_set_any(img,hdr,rule,comment);
 hd->vtype=FITS_VDOUBLE;
 hd->vdouble=val;
 return(0);
}
int fits_headerset_set_string(fitsheaderset *img,char *hdr,int rule,char *str,char *comment)
{
 fitsheader	*hd;

 hd=fits_headerset_set_any(img,hdr,rule,comment);
 hd->vtype=FITS_VSTR;

 strncpy(hd->vstr,str,69);
 hd->vstr[68]=0;

 return(0);
}
int fits_headerset_set_boolean(fitsheaderset *img,char *hdr,int rule,int vbool,char *comment)
{
 fitsheader	*hd;
 hd=fits_headerset_set_any(img,hdr,rule,comment);
 hd->vtype=FITS_VBOOLEAN;
 hd->vint=vbool;
 return(0);
}

/*****************************************************************************/

int fits_headerset_delete(fitsheaderset *img,char *hdr,int k)
{
 int	id,n;
 n=fits_headerset_get_count(img,hdr);
 if ( n<=0 || k>=n )	return(1);
 id=fits_headerset_get_id(img,hdr,k);
 memmove(img->hdrs+id,img->hdrs+id+1,(img->nhdr-id-1)*sizeof(fitsheader));
 img->nhdr--;
 return(0);
}
int fits_headerset_delete_all(fitsheaderset *img,char *hdr)
{
 int	id,n;

 n=fits_headerset_get_count(img,hdr);
 while ( n>0 )
  {	n--;
	id=fits_headerset_get_id(img,hdr,n);
	memmove(img->hdrs+id,img->hdrs+id+1,(img->nhdr-id-1)*sizeof(fitsheader));
	img->nhdr--;
  }
 return(0);
}

/*****************************************************************************/

int fits_headerset_copy(fitsheaderset *im1,fitsheaderset *im2)
{
 if ( im1->hdrs != NULL )	free(im1->hdrs);
 im1->nhdr=im1->ahdr=0;

 im1->hdrs=(fitsheader *)malloc(sizeof(fitsheader)*im2->ahdr);
 memcpy(im1->hdrs,im2->hdrs,sizeof(fitsheader)*im2->nhdr);
 im1->nhdr=im2->nhdr;
 im1->ahdr=im2->ahdr;
 return(0);
}

/*****************************************************************************/

int fits_headerset_merge(fitsheaderset *ximg,fitsheaderset *img,int inherit)
{
 int		i,id,cnt,is_inherit;
 fitsheader	*hdr;

 fits_headerset_delete_all(ximg,headers[HDR_XTENSION]);
 fits_headerset_set_boolean(ximg,headers[HDR_EXTEND],FITS_SH_INSERT,0,comment_fits_standard);
 fits_headerset_set_boolean(ximg,headers[HDR_SIMPLE],FITS_SH_INSERT,1,comment_fits_standard);

 if ( inherit<0 )	return(0);
 
 cnt=fits_headerset_get_count(ximg,headers[HDR_INHERIT]);
 if ( cnt>0 )	hdr=fits_headerset_get_header(ximg,headers[HDR_INHERIT],0);
 else		hdr=NULL;
 is_inherit=(hdr != NULL && hdr->vtype==FITS_VBOOLEAN && hdr->vint && img != NULL);
 if ( ! is_inherit && ! inherit )	return(0);

 id=fits_headerset_get_id(img,headers[HDR_EXTEND],0);
 if ( id<0 )		return(0);

 for ( i=id+1 ; i<img->nhdr ; i++ )
  {	hdr=fits_headerset_append(ximg);
	if ( strcmp(img->hdrs[i].name,headers[HDR_INHERIT]) )
	 {	memcpy(hdr,&img->hdrs[i],sizeof(fitsheader));	}
  }
 fits_headerset_delete_all(ximg,headers[HDR_INHERIT]);

 return(0);
}

/*****************************************************************************/

int fits_headerset_is_extension(fitsheaderset *header)
{
 fitsheader	*hx;
 char		*p;

 if ( header->nhdr <= 0 )	return(-1);
 else if ( header->hdrs==NULL )	return(-1);

 hx=&header->hdrs[0];

 if ( hx != NULL && hx->vtype==FITS_VSTR )
  {	p=strchr(hx->vstr,' ');
	if ( p != NULL ) *p=0;

	if ( strcmp(hx->vstr,"IMAGE")==0 )
		return(FITS_EXT_IMAGE);
	else if ( strcmp(hx->vstr,"TABLE")==0 )
		return(FITS_EXT_TABLE);
	else if ( strcmp(hx->vstr,"BINTABLE")==0 )
		return(FITS_EXT_BINTABLE);
	else
		return(-1);
  }
 else	
	return(-1);
}
 
/*****************************************************************************/
