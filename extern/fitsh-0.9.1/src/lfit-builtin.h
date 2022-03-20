/*****************************************************************************/
/* lfit-builtin.h							     */
/*****************************************************************************/

#ifndef	__LFIT_BUILTIN_H_INCLUDED
#define	__LFIT_BUILTIN_H_INCLUDED	1

/*****************************************************************************/

#include <lfit/lfit.h>

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

#define		F_ARG		48
#define		F_ATAN2		49
#define		F_HYPOT		50
#define		F_VPI		51

#define		F_EXP		52
#define		F_LOG		53
#define		F_EXP10		54
#define		F_LOG10		55

#define		F_SQR		56
#define		F_ABS		57
#define		F_SGN		58
#define		F_THETA		59
#define		F_FMOD		60
#define		F_FINT		61
#define		F_FDIV		62

#define		F_ILINEAR	24
#define		F_ILINEAR_DX	25
#define		F_ISPLINE	26
#define		F_ISPLINE_DX	27
#define		F_ICYSPLINE	28
#define		F_ICYSPLINE_DX	29

#define		F_JBESSEL	64
#define		F_YBESSEL	65

/*****************************************************************************/

#define		LFIT_F_MAX_BUILTIN		128	/* incl. nonexported */

/*****************************************************************************/

extern	short *psn_lfit_simp[];

extern	psnlfit	psnlfit_list_builtin_normal_operators[];
extern	psnlfit	psnlfit_list_builtin_relational_operators[];
extern	psnlfit psnlfit_list_builtin_elementary_functions[];
extern	psnlfit	psnlfit_list_builtin_interpolators[];
extern	psnlfit	psnlfit_list_builtin_aa_functions[];
#if defined _FI_SOURCE || defined _ASTRO_EXTEN
extern	psnlfit psnlfit_list_builtin_xfuncts[];
#endif

/*****************************************************************************/

#endif

/*****************************************************************************/
                                                       
