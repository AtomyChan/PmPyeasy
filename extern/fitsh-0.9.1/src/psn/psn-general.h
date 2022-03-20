/*****************************************************************************/
/* psn_general.h							     */
/*****************************************************************************/

#ifndef __PSN_GENERAL_H_INCLUDED
#define __PSN_GENERAL_H_INCLUDED	1

/*****************************************************************************/

#include	<psn/psn.h>

/*****************************************************************************/

#define		O_ADD		1	/* + */
#define		O_SUB		2	/* - */
#define		O_MUL		3	/* * */
#define		O_DIV		4	/* / */
#define		O_POW		5	/* ^ */
#define		O_CHS		6	/* change sign, unary variant of '-' */

#define		O_RCP		7	/* Reciprocal	   // these two	     */
#define		O_PSQ		8	/* Square (power)  // have no symbol */

#define		O_EQU		10
#define		O_NEQ		11
#define		O_LS		12
#define		O_GR		13
#define		O_LE		14
#define		O_GE		15
#define		O_AND		20
#define		O_OR		21

/* Basic functions */
#define		F_SIN		32
#define		F_COS		33
#define		F_TAN		34
#define		F_CTN		35

#define		F_SINA		36
#define		F_COSA		37
#define		F_TANA		38
#define		F_CTNA		39

#define		F_ASIN		40
#define		F_ACOS		41
#define		F_ATAN		42
#define		F_ACTN		43

#define		F_ASINA		44
#define		F_ACOSA		45
#define		F_ATANA		46
#define		F_ACTNA		47

#define		F_EXP		48
#define		F_LOG		49
#define		F_EXP10		50
#define		F_LOG10		51

#define		F_SQR		52
#define		F_ABS		53
#define		F_SGN		54
#define		F_FMOD		55
#define		F_FINT		56
#define		F_FDIV		57
#define		F_ARG		58

#define		F_VPI		59

#define		F_RND		60
#define		F_GAUSS		61

/*****************************************************************************/

extern		psnsym		psn_general_op[];
extern		psnsym		psn_general_xop[];
extern		psnsym		psn_general_fn[];
extern		psnsym		psn_general_fn_rnd[];

extern		psnprop		psn_general_prop[];

extern		psnfunct	psn_general_funct[];

/*****************************************************************************/

#endif

/*****************************************************************************/
            
