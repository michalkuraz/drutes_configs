#!/bin/bash

sed -e 's/DD/'$1'/g' ../drutes.conf/ADE/contaminant.conf.temp > ../drutes.conf/ADE/contaminant.conf

# echo $2

# sed -e 's/EE/'$2'/g' ../drutes.conf/ADE/contaminant.conf.temp1 > ../drutes.conf/ADE/contaminant.conf

cd ..

bin/drutes

cd -
