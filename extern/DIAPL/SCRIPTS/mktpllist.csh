#! /bin/tcsh -f

#========================================================#
#                                                        #
#  mktpllist.csh          version 1.0.2   2006.10.23     #
#                                                        #
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

echo "   mktpllist"
echo "   ---------"
echo " "

source diapl_setup.csh

set FWHM_FILE = "fwhm.log"
set BKG_DENOM =  2
set FWHM_DENOM = 10

cd ${WORK_DIR}

#-------------------------------
foreach im (`cat ${IMAGES}`)
  ln -v -s ${DATA_DIR}/${im:gr}.${FITS}
end

echo " "
echo "FWHMs are calculated"
echo "--------------------"
rm -vf ${FWHM_FILE}
touch ${FWHM_FILE}
foreach im (`cat ${IMAGES}`)
  ${BIN}/fwhm ${im:gr}.${FITS} >> ${FWHM_FILE}
end

echo " "
echo "RM: remove symbolic links"
echo "-------------------------"
foreach im (`cat ${IMAGES}`)
  rm -v ${im:gr}.${FITS}
end

set N = `wc ${IMAGES}`
@ N = ${N[1]}
@ N1 = ${N} / ${BKG_DENOM}
@ N2 = ${N} / ${FWHM_DENOM}

sort -k 3 -n ${FWHM_FILE} | head -${N1} \
  | sort -k 2 -n | head -${N2} \
  | awk '{print $1}' > ${TPLIMAGES}

echo ${N2} image names written to ${TPLIMAGES}
echo "DONE"

exit 0
