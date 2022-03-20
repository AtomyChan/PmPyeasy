#! /bin/tcsh -f

#========================================================#
#                                                        #
#  phot.csh               version 1.0.2   2006.10.23     #
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

set names=`cat ${IMAGES}`
set images=`echo ${names:gr}`

@ xload = ( ${NX0} + ${NX} - 1 ) / ${NX}
@ yload = ( ${NY0} + ${NY} - 1 ) / ${NY}

@ nxw = ${NX} + 2 * ${MARG}
@ nyw = ${NY} + 2 * ${MARG}

echo "MAIN LOOP m by n subframes"

set images1 = `echo ${images}`
set names = `fgrep -v ${REFIM} ${IMAGES}`
set images2 = `echo ${names:gr}`

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

      set failed = 0

      @ xl -= 1
      @ yl -= 1

      echo "offsets: "${xl} ${yl}
      ${BIN}/phot phot.par ref_${prefix}.${FITS} psf_${prefix}.bin \
                  ${photimages} ${FIELD}_${prefix}
      if ( ${status} != 0 ) then
         set failed = 1
         set pname = "phot"
         goto "CRASH"
      endif

#-----------------------------

CRASH:
      if ( ${failed} != 0 ) then
         echo "CRASH for section [${ix}, ${iy}]; program: ${pname}"
      else
         echo  "completed section [${ix}, ${iy}]"
      endif

      @ ix += 1
   end
   @ ix  = 1
   @ iy += 1
end

echo END OF THE SUBFRAME LOOP

exit 0
