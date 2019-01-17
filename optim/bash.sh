#!/bin/bash

function run_drutes {

 
  rm -rf $1

  mkdir $1

  cp -a drutemp/* $1

  cd $1

  
  k=$2

  k=`echo ${k} | sed -e 's/[eE]+*/\\*10\\^/'`
 

  k=`echo "e($k)" | bc -l`

  echo $k

  sed -e 's/!KD/'$k'/g'  drutes.conf/ADE/sorption.conf.temp > drutes.conf/ADE/sorption.conf
  
  bin/drutes > /dev/null
  
  cat out/objfnc.val | grep -v "#" > objfnc.val
  
  read val < objfnc.val
  
  echo "$2 $val" >> ../drutes.vals
  
  cd .. 
    
}

rm -f drutes.vals       
let nproc=0
while read l a
  do
    if [[  $l == "p"  ]]; then
      let nproc=nproc+1
    fi
  done < pars.in
  

let z=0
while read l a
  do
    if [[  $l == "p"  ]]; then
      let z=z+1
      if [[ $z -lt $nproc ]] ; then
        run_drutes $z $a  &
      else
        run_drutes $z $a  
      fi
    fi
  done < pars.in
 

let z=0
    
while [[ $z -lt $nproc ]]
  do
  let z=0
    for i in `seq 1 $nproc`
      do
        FILE="$i/objfnc.val"
        if [ -f $FILE ];
          then
             let z=z+1
        fi
    done    
done
