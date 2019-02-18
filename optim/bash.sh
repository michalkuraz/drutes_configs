#!/bin/bash

function run_drutes {

 
  rm -rf $1

  mkdir $1

  cp -a drutemp/* $1

  cd $1

  
  ka=$2
  kd=$3


  D=$5
  
  q=$4

  ka=`echo ${ka} | sed -e 's/[eE]+*/\\*10\\^/'`
 

  ka=`echo "e($ka)" | bc -l`

  kd=`echo ${kd} | sed -e 's/[eE]+*/\\*10\\^/'`
  kd=`echo "e($kd)" | bc -l`



  sed -e 's/!KA/'$ka'/g' -e 's/!KD/'$kd'/g'   drutes.conf/ADE/sorption.conf.temp > drutes.conf/ADE/sorption.conf

  sed -e 's/!D/'$D'/g'  drutes.conf/ADE/contaminant.conf.temp > drutes.conf/ADE/contaminant.conf
  
  sed -e 's/!q/'$q'/g' drutes.conf/ADE/ADE.conf.temp > drutes.conf/ADE/ADE.conf

  
  bin/drutes -o optim > /dev/null
  
  cat out/objfnc.val | grep -v "#" > out/objfnc.val.2
  
  read val < out/objfnc.val.2

  val=`echo ${val} | sed -e 's/[eE]+*/\\*10\\^/'`
  
  val=`echo "${val}*100000" | bc -l` 
  
  echo $val > objfnc.val
  
  read oldval < ../oldval
  
  st=`echo "$val < $oldval" | bc -l`
  
  read itcount < ../itcount
  
  let itcount=$itcount+1
    
  if [[ st -eq 1 ]]; then
     echo $val > ../oldval
     echo $itcount $val $1 $2 $3 $4 $5 >> ../converge
  fi

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
