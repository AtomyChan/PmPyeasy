#include "CCDProc.h"

/***************************************************************************
 *
 * General purpose utility functions                       
 *
 * Contents:
 *    LeftStr - returns n left-most characters from a string
 *    RightStr - returns n right-most characters from a string
 *    MidStr - returns n middle characters starting at a given position
 *             within a string
 *    BZero - clears (zeros) a string
 *    BCmp - byte-wise comparison of the first length characters
 *           of two strings
 *    UpperCase - convert a string to all uppercase
 *    GetArg - general-purpose command-line argument parser
 *    ltos - convert a long integer to a string
 *
 * Shamelessly stolen and adapted from the utility routines in the 
 * Caliban program by B. Hartung.
 *
 * R. Pogge 
 * OSU Astronomy Dept
 * 1998 April 2
 *
 ***************************************************************************/

/* #include "autolog.h" */

/***************************************************************************
 * LeftStr Function                                                 
 * Purpose: Returns n left-most characters from a given string      
 * Requires: A destination buffer, the count, and the source string 
 * Returns: With dest containing n-most left characters of source   
 */

void LeftStr(char *dest, char *source, int length)
{
  int i;

  for(i=0;i<length && i<strlen(source);i++) {
    dest[i] = source[i];
  }
  dest[i] = NUL;
}

/***************************************************************************
 * RightStr Function                                                
 * Purpose: Returns n right-most characters from a given string     
 * Requires: A destination buffer, the count, and the source string 
 * Returns: With dest containing n-most left characters of source   
 */

void RightStr(char *dest, char *source, int length)
{
  int i,j;
  int sourcelen;

  if(length>(sourcelen=strlen(source)))
    strcpy(dest, source);
  else
    {
      for(i=0,j=sourcelen-length;i<length;i++,j++)
	dest[i] = source[j];
      dest[i] = NUL;
    }
}

/***************************************************************************
 * MidStr Function                                                  
 * Purpose: Returns length characters starting at start from source 
 * Requires: A destination buffer, source string, start, and length 
 * Returns: With n middle characters from source in destination     
 */

void MidStr(char *dest, char *source, int start, int length)
{
  int i;

  for(i=0,start--;i<length && start<strlen(source);i++,start++)
    dest[i]=source[start];
}

/***************************************************************************
 * BZero Function                                                   
 * Purpose: Initialize an array with nulls                          
 * Requires: A buffer and its length                                
 * Returns: With buffer containing length nulls                     
 */

void BZero(char *b, int length)
{
  for(;length>0;length--)
    b[length-1]=NUL;
}

/***************************************************************************
 * BCmp Function
 * Purpose: Performs a byte-wise comparison of the first length characters
 * of two strings 
 * Requires: Two strings and the number of characters to check 
 * Returns: TRUE or FALSE depending on if the first length characters match 
 */

int BCmp(char *b, char *c, int length)
{
  for(;length>0;length--)
    if(b[length-1] != c[length-1])
      return(FALSE);

  return(TRUE);
}

/***************************************************************************
 * UpperCase Function                                               
 * Purpose: Converts all characters in an array to uppercase        
 * Requires: A character array                                      
 * Returns: With the string having been converted to uppercase      
 */

void UpperCase(char *str)
{
  int i;

  for(i=0;i<strlen(str);i++)
    if(str[i] >= 'a' && str[i] <= 'z')
      str[i] = str[i] - ('a' - 'A');
}

/***************************************************************************
 *
 * GetArg Function
 *
 * Purpose: General purpose command-line argument parser
 * 
 * Requires: The source string of arguments, an arg number and the
 *           destination buffer
 *
 * Returns: With the given argument stored in retarg 
 * 
 */

void GetArg(char *argstr, int argnum, char *retarg)
{
  int arglen, i=0, j=0;
  
  arglen=strlen(argstr);        /* Record the length of the overall argument */
                                /* string */

  for(i=0;argnum>1;argnum--) {  /* Seek to the beginning of the nth argument */
                                /* using the space character between tokens  */
    while(argstr[i]!=' ') {
      i++;
      if(i>=arglen)             /* If you hit the end, there aren't at least */
				/* n arguments, so return null */
	break;
    }

    while(argstr[i]==' ') {     /* Ignore extra spaces between args */
      i++;
      if(i>=arglen) {
	break;
      }
    }
  }

  /* Now we should be on the first character of the nth argument So, begin
   * copying it to the return buffer until the next space or end of the arg
   * string.  Replace any nulls with a space character. */

  for(j=0;(i<arglen) && ((argstr[i]!=' ') || 
			 (argstr[i]=='"')) && (argstr[i]!='\n'); i++, j++) {
    retarg[j] = (argstr[i]==NUL ? ' ' : argstr[i]); 
  }

  retarg[j] = NUL;
}

/***************************************************************************
 * ltos Function 
 * Purpose: Converts a long to a string 
 * Requires: A string buffer and a long 
 * Returns: With the string representation of long in argbuf 
 */

void ltos(char *argbuf, long num) {
  int i=0, j=0;
  char c;
  char buf[SHORT_STR_SIZE];

  /* The algorithm used chops one digit off at a time (via mod 10) and
   * converts it into a character (via adding '0').  This series of
   * characters gets written into buf, but in reverse order, necessitating
   * the final for loop which reverses them back */

  while(num>0) {
    c=(num%10)+'0';
    num=num/10;
    buf[i++]=c;
  }
  argbuf[i] = NUL;
  
  for(j=0;i>0;argbuf[j++]=buf[--i]);
}

