/*****************************************************************************/
/* star-model.h								     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Auxiliary definitions for star-model.c, not included by other modules.    */
/*****************************************************************************/

#ifndef	__STAR_MODEL_C
#error	"This file (star-model.h) can only be included from star-model.c"
#endif

#ifndef	__STAR_MODEL_H_INCLUDED
#define	__STAR_MODEL_H_INCLUDED	1

#define		MXO	MAX_DEVIATION_ORDER		/* =  4	*/
#define		MXC	MAX_DEVIATION_COEFF		/* = 15	*/
#define		MXT	(MXC+MXO+2+MXO+3)		/* = 28	*/

#define		I00	 0
#define		I10	 1
#define		I01	 2
#define		I20	 3
#define		I11	 4
#define		I02	 5
#define		I30	 6
#define		I21	 7
#define		I12	 8
#define		I03	 9
#define		I40	10
#define		I31	11
#define		I22	12
#define		I13	13
#define		I04	14
#define		I50	15
#define		I41	16
#define		I32	17
#define		I23	18
#define		I14	19
#define		I05	20
#define		I60	21
#define		I51	22
#define		I42	23
#define		I33	24
#define		I24	25
#define		I15	26
#define		I06	27

#define		M00	(1.0)
#define		M10	(1.0)
#define		M01	(1.0)
#define		M20	(1.0/2.0)
#define		M11	(1.0)
#define		M02	(1.0/2.0)
#define		M30	(1.0/6.0)
#define		M21	(1.0/2.0)
#define		M12	(1.0/2.0)
#define		M03	(1.0/6.0)
#define		M40	(1.0/24.0)
#define		M31	(1.0/6.0)
#define		M22	(1.0/4.0)
#define		M13	(1.0/6.0)
#define		M04	(1.0/24.0)

#endif

