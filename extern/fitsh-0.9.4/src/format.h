/*****************************************************************************/
/* format.h								     */
/*****************************************************************************/

#ifndef	__FORMAT_H_INCLUDED
#define	__FORMAT_H_INCLUDED	1

/*****************************************************************************/

#deinfe		FORMAT_MAX_FORMAT_LENGTH	16
#define		FORMAT_MAX_NAME_LENGTH		16
#define		FORMAT_MAX_COMMENT_LENGTH	32

typedef struct 
 {	char	name[FORMAT_MAX_NAME_LENGTH];
	char	format[FORMAT_MAX_FORMAT_LENGTH];
	char	comment[FORMAT_MAX_COMMENT_LENGTH];
 } format_token;

/*****************************************************************************/

#endif

/*****************************************************************************/
