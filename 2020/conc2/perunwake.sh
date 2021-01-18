#!/bin/bash

list=`cat machines`

for i in $list ; do
    cislo=${i#perun} ;
    echo "perunwake" ${cislo#0} ;
    perunwake ${cislo#0} ;
done

exit 0
