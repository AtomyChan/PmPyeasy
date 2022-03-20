#! /bin/tcsh -f

#========================================================#
#                                                        #
#  shifts.csh             version 1.6.2   2006.10.23     #
#                                                        #
#  Copyright (C) 2006 by Przemek Wozniak                 #
#  wozniak@lanl.gov                                      #
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

echo "   shifts"
echo "   ------"
echo " "

source diapl_setup.csh

cd ${WORK_DIR}

#-------------------------------
ln -v -s ${DATA_DIR}/${REFIM}.${FITS}
foreach im (`cat ${IMAGES}`)
  ln -v -s ${DATA_DIR}/${im:gr}.${FITS}
end

echo " "
echo "CROSS: shifts are calculated"
echo "-------------------------------"
${BIN}/cross cross.par ${INSTRUMENT} ${SHIFTS} ${REFIM}.${FITS} -f ${IMAGES}

echo " "
echo "RM: remove symbolic links"
echo "-------------------------------"
rm -v ${REFIM}.${FITS}
foreach im (`cat ${IMAGES}`)
  rm -v ${im:gr}.${FITS}
end

echo "DONE"

exit 0
