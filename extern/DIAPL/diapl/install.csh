#! /bin/tcsh -f
#
# Copyright 2006 by Wojtek Pych, CAMK PAN
#
# Version 1.4	2006.10.25
#
mkdir -v DIAPL_BIN
#
cd aga
pwd
make
if ($status) exit(-1)
mv -v aga ../DIAPL_BIN/
rm -vf *.o
cd ..
echo "---"
#
cd cross
pwd
make
if ($status) exit(-1)
mv -v cross ../DIAPL_BIN/
rm -vf *.o
cd ..
echo "---"
#
cd cutfitsim
pwd
make
if ($status) exit(-1)
mv -v cutfitsim ../DIAPL_BIN/
rm -vf *.o
cd ..
echo "---"
#
cd fitshedit
pwd
make
if ($status) exit(-1)
mv -v fitshedit ../DIAPL_BIN/
rm -vf *.o
cd ..
echo "---"
#
cd fwhm
pwd
make
if ($status) exit(-1)
mv -v fwhm ../DIAPL_BIN/
rm -vf *.o
cd ..
echo "---"
#
cd getpsf
pwd
make
if ($status) exit(-1)
mv -v getpsf ../DIAPL_BIN/
rm -vf *.o
cd ..
echo "---"
#
cd getvar
pwd
make
if ($status) exit(-1)
mv -v getvar ../DIAPL_BIN/
rm -vf *.o
cd ..
echo "---"
#
cd im2float
pwd
make
if ($status) exit(-1)
mv -v im2float ../DIAPL_BIN/
rm -vf *.o
cd ..
echo "---"
#
cd mkushortone
pwd
make
if ($status) exit(-1)
mv -v mkushortone ../DIAPL_BIN/
rm -vf *.o
cd ..
echo "---"
#
cd mstack
pwd
make
if ($status) exit(-1)
mv -v mstack ../DIAPL_BIN/
rm -vf *.o
cd ..
echo "---"
#
cd phot
pwd
make
if ($status) exit(-1)
mv -v phot ../DIAPL_BIN/
rm -vf *.o
cd ..
echo "---"
#
cd resample2
pwd
make
if ($status) exit(-1)
mv -v resample2 ../DIAPL_BIN/
rm -vf *.o
cd ..
echo "---"
#
cd sfind
pwd
make
if ($status) exit(-1)
mv -v sfind ../DIAPL_BIN/
rm -vf *.o
cd ..
echo "---"
#
cd template
pwd
make
if ($status) exit(-1)
mv -v template ../DIAPL_BIN/
rm -vf *.o
cd ..
echo "---"
#
cd xygrid
pwd
make
if ($status) exit(-1)
mv -v xygrid ../DIAPL_BIN/
rm -vf *.o
cd ..
echo "---"
#
cd xymatch
pwd
make
if ($status) exit(-1)
mv -v xymatch ../DIAPL_BIN/
rm -vf *.o
cd ..
echo "---"
#
echo ""
echo "        DIAPL: Compilation finished"
echo ""
