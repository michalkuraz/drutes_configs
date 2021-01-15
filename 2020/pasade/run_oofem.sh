#!/bin/bash
#
# This script will run oofem
#
# M. Leps  2007
#

if [ "$1" = "" ] ; then
    echo
    echo "Usage: run_oofem.sh input_file.in"
    echo
    exit 1
fi

temp="${1%.in}"

./oofem -f $1 > ${1%.in}.out2 2>&1 
./extractor -f $1 >> ${1%.in}.dat 2>/dev/null
perl dotime3.pl < ${1%.in}.out > ${1%.in}.time

#echo "rm ${temp}.out" "${temp}.out2" "${temp}.*.osf"
rm -f "${temp}.out" "${temp}.out2" *.osf

exit 0