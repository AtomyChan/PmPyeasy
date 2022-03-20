#========================================================#
#                                                        #
#  diapl_setup.csh        version 1.7.0   2005.10.27     #
#                                                        #
#  Copyright (C) 2005 by Wojtek Pych, CAMK PAN           #
#  pych@camk.edu.pl                                      #
#                                                        #
# To be source'd by:                                     #
#   shifts.csh, template.csh and pipe.csh                #
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

#--------------------------------------------------------#
# DATA_DIR:   Original images directory
# WORK_DIR:   Working directory: scripts, parameter files, results
# BIN:        Executable binaries directory
# INSTRUMENT: Instrument parameters file
# IMAGES:     List of the images
# TPLIMAGES:  List of template images
# SHIFTS:     Shifts file
# REFIM:      Coordinates reference image
# TPLNAME:    Prefix of the filename for the template image.
# EDGE:       "Limit for shifts versus transformation"
# FITS:       Default filename extension for FITS files.
#--------------------------------------------------------#

set DATA_DIR   = ""
set WORK_DIR   = ""
set BIN        = "DIAPL_BIN"
set INSTRUMENT = "instrument.par"
set IMAGES     = "images.txt"
set TPLIMAGES  = "tplimages.txt"
set SHIFTS     = "shifts.txt"
set REFIM      = "im0001.fits"
set TPLNAME    = "ref"
set EDGE       = 5

#--------------------------------------------------------#
# Do not edit below this line                            #
#--------------------------------------------------------#

set REFIM = ${REFIM:r}

set FITS = `grep -w FITS ${INSTRUMENT}`
set FITS = ${FITS[3]}

# Input image size
set NX0 = `${BIN}/fitshedit ${DATA_DIR}/${REFIM}.${FITS} NAXIS1`
set NY0 = `${BIN}/fitshedit ${DATA_DIR}/${REFIM}.${FITS} NAXIS2`

set NX0 = ${NX0[3]}
set NY0 = ${NY0[3]}

# Sub-Image size 
set NX = `grep -w NX ${INSTRUMENT}`
set NY = `grep -w NY ${INSTRUMENT}`
set MARG = `grep -w MARG ${INSTRUMENT}`

set NX = ${NX[3]}
set NY = ${NY[3]}
set MARG = ${MARG[3]}
