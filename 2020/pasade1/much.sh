#!/bin/sh
#
# This script will run `command [params]`
# "number" times.
#
# O. Kupca 2001
# M. Leps  2002, 2004
#

if [ "$1" = "" ] || [ "$2" = "" ]; then
    echo
    echo "Usage: much.sh number command [params] "
    echo
    exit 1
fi

count=0
first=$1
shift
name=$1
shift

while [ $count -lt $first ]
do
    count=$(($count  + 1))
    echo $count
    $name $count $@
done

echo "*********************"
echo "***** *    * ****"
echo "*     **   * *   *"
echo "***   * *  * *   *"
echo "*     *  * * *   *"
echo "*     *   ** *   *"
echo "****  *    * ****"
echo "********************"

exit 0
