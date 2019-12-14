#!/bin/bash

let j=0
for i in `ls inputs/go | grep dat` ; do
  let j=$j+1
  cp inputs/go/$i drutemp/drutes.conf/kinwave/inputs.dat
  echo "processing datafiles $i"

  while read ch val ; do
    if [[  $ch == "i"  ]]; then
      rain=$val
    fi
    if [[ $ch == "s" ]] ; then
      slope=$val
    fi
  done < inputs/go/$i.in

  sed -e 's/!tan/'$slope'/g' drutemp/drutes.conf/kinwave/kinwave.conf.init > drutemp/drutes.conf/kinwave/kinwave.conf.temp
  sed -e 's/!R/'$rain'/g' drutemp/drutes.conf/kinwave/rain.init > drutemp/drutes.conf/kinwave/rain.in
  
  rm -rf results

  rm -rf extremes.sav 

  rm -rf totvals


  echo "running SADE optim for set $i"  

  ./pasade 1>all_log.out 2>&1
  

  echo "SADE has finished"

  rm -rf pars.in

  rm -rf counter
 
  if [  -s extremes.sav ] ; then

    while read count S A N ; do
      echo "t" $S $A $N  >> pars.in
      echo $count > counter
    done < extremes.sav

  fi

  read null1 null2 null3 null4 S A N < results
  echo "t" $S $A $N  >> pars.in

  if [ -f counter ]; then
    read count < counter
    let count=$count+1
  else
    let count=1
  fi


  echo "number of extremes was: $count"

  ./bash.sh
 
 
  if [ ! -d outputs ] ; then 
    mkdir outputs
  fi
  
  echo "packing all outputs..."
  tar -czf outputs/$i.out.tgz $(seq $count) totvals extremes.sav results
 
done
