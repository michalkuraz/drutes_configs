#!/bin/bash

function run_drutes {

  rm -rf $1

  mkdir $1

  cp -a drutemp/* $1

  cd $1
  
  a=$2
  n=$3
  ths=$4
  Ks=$5
  Ss=$6

  m=`echo "scale=12; 1.0-1.0/$n" | bc`
   
  sed -e 's/!A/'$a'/g' -e 's/!N/'$n'/g' -e 's/!M/'$m'/g' -e 's/!K/'$Ks'/g' -e 's/!T/'$ths'/g' -e 's/!S/'$Ss'/g' drutes.conf/water.conf/matrix.conf.temp > drutes.conf/water.conf/matrix.conf
  
  bin/drutes -o optim > /dev/null

  cp out/objfnc.val ./
  
#  Rscript ../evol.R
  
 # read val < objfnc.val

  cat objfnc.val >> ../drutes.vals

  cd ..
    
}

rm -f drutes.vals       
let nproc=0
while read l a b c d e
  do
    if [[  $l == "p"  ]]; then
      let nproc=nproc+1
    fi
  done < pars.in
  

let z=0
while read l a b c d e
  do
    if [[  $l == "p"  ]]; then
      let z=z+1
      if [[ $z -lt $nproc ]] ; then
        run_drutes $z $a $b $c $d $e &
      else
        run_drutes $z $a $b $c $d $e 
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
  
cat drutes.vals >> all.vals  
echo "----------------------" >> all.vals
