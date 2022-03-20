#! /bin/tcsh -f

#========================================================#
#                                                        #
#  template.csh           version 2.3.2    2006.10.23    #
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

echo "   template"
echo "   --------"

source diapl_setup.csh

set TRIMFILES   = "trimfiles.txt"
set STACKIMAGES = "stackimages.txt"

cd ${WORK_DIR}
pwd

rm -fv ${TRIMFILES} ${STACKIMAGES}
touch ${TRIMFILES} ${STACKIMAGES}

set NAMES=`awk '{print $1}' ${TPLIMAGES}`
if (${status} != 0) exit 1
foreach name ( ${NAMES:gr} )
  echo ${name}r.${FITS} >> ${STACKIMAGES}
end

set NAMES=`awk '{print $1}' ${TPLIMAGES} | grep -v ${REFIM}`
if (${status} != 0) exit 2
set IMAGES = `echo ${NAMES:gr}`

@ xload = ( ${NX0} + ${NX} - 1 ) / ${NX}
@ yload = ( ${NY0} + ${NY} - 1 ) / ${NY}

@ nxw = ${NX} + 2 * ${MARG}
@ nyw = ${NY} + 2 * ${MARG}

# SET SECTION FROM WHICH WE START
set ix = 1
set iy = 1

echo ""
echo "MAIN LOOP "${xload}" by "${yload}" subframes"
echo "=========================="

while( ${iy} <= ${yload} )
  while( ${ix} <= ${xload} )

    set prefix = "${ix}_${iy}"

    @ xl = 1 + ${NX} * (${ix} - 1) - ${MARG}
    @ xu =     ${NX} *  ${ix}      + ${MARG}
    @ yl = 1 + ${NY} * (${iy} - 1) - ${MARG}
    @ yu =     ${NY} *  ${iy}      + ${MARG}

    echo ""
    echo "--------------------------"
    echo "Section:   "${prefix}"   [ "${xl} : ${xu} , ${yl} : ${yu} ]
    echo "--------------------------"
    echo ""
    echo "FIND STARS ON THE IMAGES (sfind)"
    echo ""

# *** Reference image ***

    set arg = ` grep ${REFIM} ${SHIFTS} `
    if (${status} != 0) then
      echo "ERROR: grep "${REFIM} ${SHIFTS}": failed"
      exit 3
    endif

    @ n11 = ${xl} - ${arg[2]}
    @ n12 = ${xu} - ${arg[2]}
    @ n21 = ${yl} - ${arg[3]}
    @ n22 = ${yu} - ${arg[3]}

    echo ${REFIM}: ${n11} ${n12} ${n21} ${n22}

    ${BIN}/cutfitsim ${DATA_DIR}/${REFIM}.${FITS} ${REFIM}c.${FITS} \
                     ${n11} ${n12} ${n21} ${n22}
    if (${status} != 0) then
      echo "cutfitsim failed - aborting now"
      exit 4
    endif

    ${BIN}/sfind sfind.par ${INSTRUMENT} ${REFIM}c.${FITS} ${REFIM}c.coo
    if (${status} != 0) then
       echo "sfind failed - aborting now"
       exit 5
    endif

    ${BIN}/im2float ${REFIM}c.${FITS} ${REFIM}r.${FITS}
# *** Reference image ***

    foreach im ( ${IMAGES} )
      set arg = ` fgrep ${im} ${SHIFTS} `
      if (${status} != 0) then
        echo "ERROR: fgrep "${im} ${SHIFTS}": failed"
        exit 3
      endif

      @ n11 = ${xl} - ${arg[2]}
      @ n12 = ${xu} - ${arg[2]}
      @ n21 = ${yl} - ${arg[3]}
      @ n22 = ${yu} - ${arg[3]}

      echo ${im}: ${n11} ${n12} ${n21} ${n22}

      ${BIN}/cutfitsim ${DATA_DIR}/${im}.${FITS} ${im}c.${FITS} \
                       ${n11} ${n12} ${n21} ${n22}
      if (${status} != 0) then
        echo "cutfitsim failed - aborting now"
        exit 4
      endif

      ${BIN}/sfind sfind.par ${INSTRUMENT} ${im}c.${FITS} ${im}c.coo
      if (${status} != 0) then
         echo "sfind failed - aborting now"
         exit 5
      endif

    end

    echo ""
    echo "MATCH COORDINATES OF THE STARS (xymatch)"
    echo ""

    foreach im ( ${IMAGES} )
      echo ${im}
      set failed = 0
      set nstars = ` wc ${im}c.coo `

      if( ${nstars[1]} < 10  ||  ${status} != 0 ) then
        set failed = 1
        set pname = "sfind or wc"
        goto "FAKE"
      endif

      ${BIN}/xymatch xymatch.par ${REFIM}c.coo ${im}c.coo ${im}c.match corr.txt

      if( ${status} != 0 ) then
        set failed = 1
        set pname = "xymatch"
        goto "FAKE"
      endif

      set corr=`cat corr.txt`
      echo "xymatch -> shifts: "${corr}
      rm -f corr.txt

      if( ${corr[1]}>${EDGE} || ${corr[2]}>${EDGE} || \
          ${corr[1]}<-${EDGE} || ${corr[2]}<-${EDGE} ) then

        echo "COORECT - SHIFT GREATER THAN VALUE OF edge"

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
          exit 6
        endif

        echo ${BIN}/sfind sfind.par ${INSTRUMENT} ${im}c.${FITS} ${im}c.coo
        ${BIN}/sfind sfind.par ${INSTRUMENT} ${im}c.${FITS} ${im}c.coo
        if (${status} != 0) then
          echo "sfind failed - aborting now"
          exit 7
        endif

        set nstars = ` wc ${im}c.coo `

        if( ${nstars[1]} < 10  ||  ${status} != 0 ) then
          set failed = 1
          set pname = "sfind or wc"
          goto "FAKE"
        endif

        ${BIN}/xymatch xymatch.par ${REFIM}c.coo ${im}c.coo ${im}c.match \
                       corr.txt

        if( ${status} != 0 ) then
          set failed = 1
          set pname = "xymatch"
          goto "FAKE"
        endif
        rm -r corr.txt

      endif

      echo ""
      echo "PIXEL GRID TRANSFORMATION COEFFICIENTS ARE CALCULATED (xygrid)"
      echo ""

      ${BIN}/xygrid xygrid.par ${im}c.match ${im}c.coeff

      if( ${status} != 0 ) then
        set failed = 1
        set pname = "xygrid"
        goto "FAKE"
      endif

      echo ""
      echo "THE PIXEL GRID IS RESAMPLED TO MATCH refim PIXEL GRID (resample2)"
      echo ""

      ${BIN}/resample2 resample2.par ${INSTRUMENT} \
        ${im}c.coeff ${im}c.${FITS} ${im}r.${FITS}

      if( ${status} != 0 ) then
        set failed = 1
        set pname = "resample"
        goto "FAKE"
      endif

#-------------------------------------------------

FAKE:

      @ n11 = -${nxw}
      @ n12 = -1
      @ n21 = -${nyw}
      @ n22 = -1

      if( ${failed} != 0 ) then

        echo "FAILURE"
        echo "SET COORDINATES OUTSIDE IMAGE"

        echo  ${pname} failed, preparing fake image ${im}r.${FITS}
        ${BIN}/cutfitsim ${REFIM}r.${FITS}  ${im}r.${FITS} \
                         ${n11} ${n12} ${n21} ${n22}
      endif

#-------------------------------------------------

    end

    echo ""
    echo "END OF LOOP FOR TEMPLATE IMAGES"
    echo ""
    echo "COADDING OF THE FRAMES (mstack)"
    echo ""

    ${BIN}/mstack mstack.par ${INSTRUMENT} ${REFIM}r.${FITS} \
                  ${STACKIMAGES} ${TPLNAME}_${prefix}.${FITS}
    if (${status} != 0) then
      echo "mstack failed - aborting now"
      exit 8
    endif

    @ xl = 1     + ${MARG}
    @ xu = ${NX} + ${MARG}
    @ yl = 1     + ${MARG}
    @ yu = ${NY} + ${MARG}

    echo ""
    echo "TRIM SUBFRAMES "${xl} ${xu} ${yl} ${yu}
    echo ""

    ${BIN}/cutfitsim ${TPLNAME}_${prefix}.${FITS} trim_${prefix}.${FITS} \
                     ${xl} ${xu} ${yl} ${yu}
    if (${status} != 0) then
      echo "cutfitsim failed - aborting now"
      exit 9
    endif

    echo trim_${prefix}.${FITS} >> ${TRIMFILES}

    echo "CALCULATE THE PSF COEFFICIENTS ON A GIVEN SUBFRAME (getpsf)"
    echo ""

    ${BIN}/getpsf getpsf.par ${INSTRUMENT} ${TPLNAME}_${prefix}.${FITS} \
      psf_${prefix}.bin

    if( ${status} != 0 ) then
      echo "\n\tgetpsf failed for ${TPLNAME}_"${prefix}"."${FITS}"\n"
    endif

    @ ix += 1

  end

  @ ix  = 1
  @ iy += 1

end

echo "END OF LOOP FOR SUBFRAMES"
echo ""
echo "PREPARING A TEMPLATE IMAGE (template)"

${BIN}/template ${TRIMFILES} tpl.${FITS} ${xload} ${yload}
if (${status} != 0) then
   echo "template failed status= "${status}
endif

echo "Removing unnecessary files"

foreach f (`cat ${TRIMFILES}`)
  rm -v ${f}
end
rm -v ${TRIMFILES}

foreach f (`cat ${STACKIMAGES}`)
  rm -v ${f}
end
rm -v ${STACKIMAGES}

rm -v *.coo *.match *.coeff
rm -v *c.${FITS}

echo "DONE"

exit 0
