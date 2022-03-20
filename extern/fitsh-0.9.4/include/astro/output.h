/*****************************************************************************/
/* output.h 								     */
/*****************************************************************************/

#ifndef	__OUTPUT_H_INCLUDED
#define	__OUTPUT_H_INCLUDED	1

#define		STRDEG_SIGN		0x01
#define		STRDEG_SPACE		0x00
#define		STRDEG_COLON		0x02

/* strdeg():
   (S)prints the angle 'a' into the buffer 'buff'. The formatting parameters
   are specified by 'dl', 'prec' and 'type'. Some examples (a=1.65432):
	dl	prec	output
	0	5	1.65432
	2	2	1 39 15.55
	1	0	1 39
	1	1	1 39.2
   If 'type'&STRDEG_SIGN is non-zero the sign '+' is forced to be written
   for positive angles aslo. If 'type'&STRDEG_COLON is non-zero, the fields
   are separated with ':' instead of spaces.				     */
char *	strdeg(char *buff,double a,int dl,int prec,int type);

#endif
                                                               
