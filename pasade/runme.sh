#!/bin/bash

let j=0
for i in `ls inputs/go | grep dat` ; do
  let j=$j+1
  cp inputs/go/$i drutemp/drutes.conf/kinwave/inputs.dat
  echo "processing datafiles $j"
  let setrain=0
  while read ch val ; do
    if [[  $setrain == 0  ]]; then
      let setrain=1
      rain=$val
    else
      slope=$val
    fi
  done < inputs/go/$j.in
  
  sed -e 's/!tan/'$slope'/g' drutemp/drutes.conf/kinwave/kinwave.conf.init > drutemp/drutes.conf/kinwave/kinwave.conf.temp
  sed -e 's/!R/'$rain'/g' drutemp/drutes.conf/kinwave/rain.init > drutemp/drutes.conf/kinwave/rain.in
  
  rm -rf results
  
  rm -rf totvals
  
  ./pasade 1>all_log.out 2>&1
  
  rm -rf pars.in
  
  let count=0
  while read null1 null2 null3 null4 S A N ; do
    echo "t" $S $A $N >> pars.in
    let count=$count+1
  done < results
  
  ./bash.sh
  
  tar -czf $(seq $count) outputs/$i.out.tgz 

done
