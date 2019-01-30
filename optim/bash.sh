#!/bin/bash

function run_drutes {

 
  rm -rf $1

  mkdir $1

  cp -a drutemp/* $1

  cd $1

  
  ka=$2
  kd=$3
  cm=$4

  D=$5

  ka=`echo ${ka} | sed -e 's/[eE]+*/\\*10\\^/'`
 

  ka=`echo "e($ka)" | bc -l`

  kd=`echo ${kd} | sed -e 's/[eE]+*/\\*10\\^/'`
  kd=`echo "e($kd)" | bc -l`

  cm=`echo ${cm} | sed -e 's/[eE]+*/\\*10\\^/'`
  cm=`echo "e($cm)" | bc -l`


  sed -e 's/!KA/'$ka'/g' -e 's/!KD/'$kd'/g' -e 's/!CM/'$cm'/g'  drutes.conf/ADE/sorption.conf.temp > drutes.conf/ADE/sorption.conf

  sed -e 's/!D/'$D'/g'  drutes.conf/ADE/contaminant.conf.temp > drutes.conf/ADE/contaminant.conf

  
  bin/drutes -o optim > /dev/null
  
  cat out/objfnc.val | grep -v "#" > objfnc.val
  
  read val < objfnc.val

  val=`echo ${val} | sed -e 's/[eE]+*/\\*10\\^/'`
  
  val=`echo "${val}*100000" | bc -l` 
  
  echo $val > objfnc.val

  cd .. 
    
}

      
let nproc=0

while read l a b c d
  do
    if [[  $l == "p"  ]]; then
      let nproc=nproc+1
    fi
  done < pars.in
  


let z=0
while read l a b c d
  do
    if [[  $l == "p"  ]]; then
      let z=z+1
      if [[ $z -lt $nproc ]] ; then
        run_drutes $z $a $b $c $d  &
      else
        run_drutes $z $a $b $c $d
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
