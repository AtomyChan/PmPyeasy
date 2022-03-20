/*****************************************************************************/
/* astro.h								     */
/*****************************************************************************/

#ifndef	__ASTRO_H_INCLUDED
#define	__ASTRO_H_INCLUDED	1

/* global constants */ /******************************************************/

extern	const	double	au_km;				/* AU/km	     */
extern	const	double	r_earth_pol,r_earth_equ;	/* polar & ekvatorial*/
extern	const	double	r_sun,r_moon;			/* in kilometers...  */

extern	const	double	m_mercury,		/* mass of the planets	*/
			m_venus,		/* in solar mass units	*/
			m_earth,
			m_mars,
			m_jupiter,		
			m_saturn,
			m_uranus,
			m_neptune;

/* other header files */ /****************************************************/

#include "cmath.h"
#include "dtc.h"
#include "earth.h"
#include "easpec.h"
#include "planets.h"

/*****************************************************************************/

#endif

/*****************************************************************************/
                                         
                              

