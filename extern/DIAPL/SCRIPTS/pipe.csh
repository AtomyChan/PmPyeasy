#! /bin/tcsh -f

#========================================================#
#                                                        #
#  pipe.csh               version 1.9.5   2006.10.23     #
#                                                        #
#  Original source code:                                 #
#  Copyright (C) 2006 by Przemek Wozniak                 #
#  wozniak@lanl.gov                                      #
#                                                        #
#  Modifications:                                        #
#  Copyright (C) 2006 by Wojtek Pych, CAMK PAN           #
#  pych@camk.edu.pl                                      #
#                                                        #
#========================================================#

#*************************************************************************#
#                                                                         #
#   This program is free software; you can redistribute it and/or modify  #
#   it under the terms of the GNU General Public License as published by  #
#   the Free Software Foundation; either version 2 of the License, or     #
#   (at your option) any later version.                                   #
#                                                                         #
#   This program is distributed in the hope that it will be useful,       #
#   but WITHOUT ANY WARRANTY; without even the implied warranty of        #
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         #
#   GNU General Public License for more details.                          #
#                                                                         #
#   You should have received a copy of the GNU General Public License     #
#   along with this program; if not, write to the                         #
#   Free Software Foundation, Inc.,                                       #
#   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             #
#                                                                         #
#*************************************************************************#

source diapl_setup.csh

set FIELD="default"

set mrjimages="mrjimages.txt"
set difimages="difimages.txt"
set photimages="photimages.txt"

#-----------------------------

cd ${WORK_DIR}

echo -n "PREPARATION OF THE LISTS WITH IMAGE NAMES AND RESETTING"
echo " OF THE DATABASE CONTENTS (IF EXIST)"

set names=`cat ${IMAGES}`
set images=`echo ${names:gr}`

@ xload = ( ${NX0} + ${NX} - 1 ) / ${NX}
@ yload = ( ${NY0} + ${NY} - 1 ) / ${NY}

@ nxw = ${NX} + 2 * ${MARG}
@ nyw = ${NY} + 2 * ${MARG}

echo "SET MASK OF BAD PIXELS - HERE NOT DEFINED"

@ nmax = ${nxw}
if (${nmax} < ${nyw}) then 
  @ nmax = ${nyw}
endif

${BIN}/mkushortone ${nmax} ONE_ushort.${FITS}

${BIN}/cutfitsim ONE_ushort.${FITS} blank_ushort.${FITS} 1 ${nxw} 1 ${nyw}
if (${status} != 0) then
   echo "cutfitsim failed - aborting now"
   exit 1
endif

echo "MAIN LOOP "${xload}" by "${yload}" subframes"

set images1 = `echo ${images}`
set names = `cat ${IMAGES} | fgrep -v ${REFIM}`
set images2 = `echo ${names:gr}`

echo "SET SECTION FROM WHICH WE START"

set ix = 1
set iy = 1

while( ${iy} <= ${yload} )
   while( ${ix} <= ${xload} )

      set prefix = "${ix}_${iy}"
      echo "------------"
      echo "    "${prefix}
      echo "------------"

      @ xl = 1 + ${NX} * (${ix} - 1) - ${MARG}
      @ xu =     ${NX} *  ${ix}      + ${MARG}
      @ yl = 1 + ${NY} * (${iy} - 1) - ${MARG}
      @ yu =     ${NY} *  ${iy}      + ${MARG}


      rm -vf ${mrjimages} ${difimages} ${photimages}
      touch ${mrjimages} ${difimages} ${photimages}

      foreach name (`echo ${images}`)
        echo ${name}r${prefix}.${FITS} >> ${mrjimages}
        echo "s_"${name}r${prefix}.${FITS} >> ${difimages}
        echo "s_"${name}r${prefix}.${FITS} ${name}r${prefix}.${FITS} \
          ${name}r${prefix}.ker >> ${photimages}
      end

      echo "SEARCH FOR STARS ON THE IMAGE"

      foreach im ( ${images1} )

         set arg = ` fgrep ${im} ${SHIFTS} `

         @ n11 = ${xl} - ${arg[2]}
         @ n12 = ${xu} - ${arg[2]}
         @ n21 = ${yl} - ${arg[3]}
         @ n22 = ${yu} - ${arg[3]}

         echo ${im}: ${n11} ${n12} ${n21} ${n22}

         ${BIN}/cutfitsim ${DATA_DIR}/${im}.${FITS} ${im}c.${FITS} \
                          ${n11} ${n12} ${n21} ${n22}
         if (${status} != 0) then
            echo "cutfitsim failed - aborting now"
            exit 1
         endif

         ${BIN}/sfind sfind.par ${INSTRUMENT} ${im}c.${FITS} ${im}c.coo
         if (${status} != 0) then
            echo "sfind failed - aborting now"
            exit 2
         endif
      end

echo  " "
echo  ${BIN}/im2float ${REFIM}c.${FITS} ${REFIM}r${prefix}.${FITS}
      ${BIN}/im2float ${REFIM}c.${FITS} ${REFIM}r${prefix}.${FITS}
echo  " "

      echo "MATCH COORDINATES OF THE STARS"
      foreach im ( ${images2} )
         echo ${im}
         set failed = 0
         set nstars = ` wc ${im}c.coo `

         if ((${nstars[1]} < 10) || (${status} != 0)) then
            echo "failed: nstars=  "${nstars[1]}
            set failed = 1
            set pname = "sfind or wc"
            goto "FAKE"
         endif

      ${BIN}/xymatch xymatch.par ${REFIM}c.coo ${im}c.coo ${im}c.match corr.txt
         if( ${status} != 0 ) then
            echo "xymatch failed"
            set failed = 1
            set pname = "xymatch"
            goto "FAKE"
         endif
         set corr=`cat corr.txt`
         rm -f corr.txt
         echo "xymatch -> shifts: "${corr}

         if ( ${corr[1]} > ${EDGE} || ${corr[2]} > ${EDGE} || \
              ${corr[1]} < -${EDGE} || ${corr[2]} < -${EDGE} ) then

            echo "CORRECT WHEN SHIFT GREATER THAN edge"

            set arg = ` fgrep ${im} ${SHIFTS} `

            @ n11 = ${xl} - ${arg[2]} + ${corr[1]}
            @ n12 = ${xu} - ${arg[2]} + ${corr[1]}
            @ n21 = ${yl} - ${arg[3]} + ${corr[2]}
            @ n22 = ${yu} - ${arg[3]} + ${corr[2]}

            echo "nn: "${n11} ${n12} ${n21} ${n22}

            ${BIN}/cutfitsim ${DATA_DIR}/${im}.${FITS} ${im}c.${FITS} \
                             ${n11} ${n12} ${n21} ${n22}
            if (${status} != 0) then
               echo "cutfitsim failed - aborting now"
               exit 3
            endif

            ${BIN}/sfind sfind.par ${INSTRUMENT} ${im}c.${FITS} ${im}c.coo

            set nstars = ` wc ${im}c.coo `

            if ( ${nstars[1]} < 10  ||  ${status} != 0 ) then
               set failed = 1
               set pname = "sfind or wc"
               goto "FAKE"
            endif

            echo "correction:"
      ${BIN}/xymatch xymatch.par ${REFIM}c.coo ${im}c.coo ${im}c.match corr.txt 
            if ( ${status} != 0 ) then
               set failed = 1
               set pname = "xymatch"
               goto "FAKE"
            endif
         rm -f corr.txt
         endif

         echo "PIXEL GRID TRANSFORMATION COEFFICIENTS ARE CALCULATED"

         ${BIN}/xygrid xygrid.par ${im}c.match ${im}c.coeff
         if ( ${status} != 0 ) then
            set failed = 1
            set pname = "xygrid"
            goto "FAKE"
         endif

         echo "THE PIXEL GRID IS RESAMPLED TO MATCH refim PIXEL GRID"

         ${BIN}/resample2 resample2.par ${INSTRUMENT} ${im}c.coeff \
                          ${im}c.${FITS} ${im}r${prefix}.${FITS}
         if( ${status} != 0 ) then
            set failed = 1
            set pname = "resample2"
            goto "FAKE"
         endif

#-----------------------------
FAKE:
         if ( ${failed} != 0 ) then
            echo "IN CASE OF FAILURE - SET COORDINATES OUTSIDE IMAGE"

            @ n11 = -${nxw}
            @ n12 = -1
            @ n21 = -${nyw}
            @ n22 = -1

            echo ${pname}" failed, preparing fake image "${im}r${prefix}.${FITS}
            ${BIN}/cutfitsim ${REFIM}r${prefix}.${FITS} \
               ${im}r${prefix}.${FITS} ${n11} ${n12} ${n21} ${n22}
         endif
      end

#-----------------------------

      echo "END OF THE LOOP FOR FRAMES"

      echo "NOW THE IMAGES ARE SUBTRACTED"
      set failed = 0

      ${BIN}/aga aga.par ${INSTRUMENT} blank_ushort.${FITS} \
		${TPLNAME}_${prefix}.${FITS} ${mrjimages}
      if ( ${status} != 0 ) then
         set failed = 1
         set pname = "aga"
         goto "CRASH"
      endif

      @ xl -= 1
      @ yl -= 1

      echo "SEARCH FOR VARIABLE OBJECTS"

      echo "offsets: "${xl} ${yl}
      ${BIN}/getvar getvar.par ${INSTRUMENT} ${TPLNAME}_${prefix}.${FITS} \
             psf_${prefix}.bin ${photimages} ${FIELD}_${prefix} ${xl} ${yl}
      if ( ${status} != 0 ) then
         set failed = 1
         set pname = "getvar"
         goto "CRASH"
      endif

      echo "FIND PHOTOMETRY OF THE VARIABLES"
      ${BIN}/phot phot.par ${INSTRUMENT} ${TPLNAME}_${prefix}.${FITS} \
                  psf_${prefix}.bin ${photimages} ${FIELD}_${prefix}
      if ( ${status} != 0 ) then
         set failed = 1
         set pname = "phot"
         goto "CRASH"
      endif

#-----------------------------

CRASH:
      if ( ${failed} != 0 ) then
         echo "WHEN SUBSTRACTION FAILS"
         echo "CRASH for section [${ix}, ${iy}]; program: ${pname}"
      endif

      if ( ${status} == 0 ) then
         echo  "completed section [${ix}, ${iy}]"
      endif

      @ ix += 1
   end
   @ ix  = 1
   @ iy += 1
end

echo END OF THE SUBFRAME LOOP

exit 0
