/*****************************************************************************/
/* dtc.h 								     */
/*****************************************************************************/

#ifndef	__DTC_H_INCLUDED
#define	__DTC_H_INCLUDED	1

/* read_julian_date():
   Reads the current Julian Day using the system function clock().	     */
double  read_julian_date(void);

/* get_julian_date(): 
   Calculates the Julian Day for the day 'ye', 'mo', 'da' (UT, at midnight). */
double  get_julian_date(int ye,int mo,int da);
/* get_calendar_date():
   Converts the Julian Day 'jd' into calendar date and returns it in the
   integers 'ye', 'mo' and the double 'da'.                                  */
void    get_calendar_date(double jd,int *ye,int *mo,double *da);

/* get_easter_day():
   Returns the day of Easter Sunday in the given year 'year' from 1st of 
   March. If the value returned is greater than 31, Easter Sunday is in April
   in the given year...                                                      */
int     get_easter_day(int year);

#endif
                                         
