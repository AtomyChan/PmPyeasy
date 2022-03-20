/*****************************************************************************/
/* str.h								     */
/*****************************************************************************/

#ifndef	__STR_H_INCLUDED
#define	__STR_H_INCLUDED	1

/*****************************************************************************/

/* strkcpy():
   It works like strncpy() but ensures that 'out' is zero-terminated 
   (therefore, the length of 'out' is always _less_ than 'max').	     */
char *	strkcpy(char *out,char *in,int max);

/* strappend():
   Appends the string 'cat' to the dynamically allocated string 'string'.  
   'string' can be NULL, in this case, the function works like strdup().     */
int	strappend(char **string,char *cat);

/* strappendf():
   It is a combination of realloc(), strcat() and sprintf(). Appends the 
   printf-formatted argument 'format' to the dynamically allocated string
   'string'. It is a bit primitive, therefore use only if long-long strings
   could not occur (otherwise it also works, but can be slow).		     */
int	strappendf(char **string,char *format,...);

/* vstrappendf(): */
int	vstrappendf(char **string,char *format,va_list ap);

/*****************************************************************************/

#endif

/*****************************************************************************/
                                                        
