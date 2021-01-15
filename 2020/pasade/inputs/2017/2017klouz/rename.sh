#!/bin/bash
# Example of converting - to _ in all files in current directory where this script is.

for i in `ls *_*`
do
 NEW=`echo $i|tr '_' '-'`
 mv $i $NEW

done
