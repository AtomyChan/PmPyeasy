#/bin/bash
#
# calling sequence
# ./rc3.sh apikey > ./script
# bash ./script

for file in $(ls ~/test_temp/*fits)
do
    cmd="python ../client.py --server http://supernova.astrometry.net/api/ --apikey $1 --upload $file"
    echo $cmd
done
