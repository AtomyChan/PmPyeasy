/*****************************************************************************/
/* floodfill.c 								     */
/*****************************************************************************/

#ifndef	__FLOODFILL_H_INCLUDED
#define	__FLOODFILL_H_INCLUDED	1

/* floodfill():
   Fills an area. The start point of the flood-filling is (x0,y0). The 
   functions getp() and putp() should read and set the pixel value 
   located at (x,y). The function getp() should return 0, if the point (x,y)
   is empty, otherwise it should return a nonzero value. The putp() function
   should fill the pixel located at (x,y) and store the information somehow 
   that this pixel has been filled (so, invoked with this (x,y) point, getp()
   should return a nonzero value). If getp(...,x0,y0) was nonzero, the function
   floodfill() wouldn't do anything. The optional argument 'param' is always
   passed to the functions getp() and putp() at all invocations.	     */

int floodfill(	int x0,int y0,
		int  (*getp)(void *param,int x,int y),
		void (*putp)(void *param,int x,int y),
		void *param);

#endif
          
