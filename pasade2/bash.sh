#!/bin/bash

function run_drutes {

  rm -rf $1

  mkdir $1

  cp -a drutemp/* $1

  cd $1
	
  A=$2
  D=$3
 
  
  tisk=$4
  

 
   
#   substitution of parameters into input files for drutes 
  sed -e 's/!LA/'$A'/g' -e 's/!LD/'$D'/g'  drutes.conf/kinwave/kinwave.conf.temp > drutes.conf/kinwave/kinwave.conf


  
  #execute drutes, send the output into "black hole"
  bin/drutes > /dev/null
  
  let val=0
  
  cat out/objfnc.val | grep -v "#" > out/objfnc.val.2
  
  read val < out/objfnc.val.2

  val=`echo ${val} | sed -e 's/[eE]+*/\\*10\\^/'`
  
  val=`echo "${val}*100000" | bc -l` 
  
  echo $val > objfnc.val
  
  
  echo $A $D $val >> ../totvals
  
  if [[ $tisk == "t" ]] ; then
    gnuplot < ../plot.gnuplot
  fi
  
  cd ..
    
}


rm -f drutes.vals


#count the number of processes       
let nproc=0
while read l a b ; do
  if [[  $l == "p"  ]]  || [[ $l == "t" ]] ; then
    let nproc=nproc+1
  fi
done < pars.in

  
#execute drutes function in parallel
let z=0
while read l a b 
  do
    if [[  $l == "p"  ]]; then
      let z=$z+1
      if [[ $z -lt $nproc ]] ; then
        run_drutes $z $a $b    &
      else
        run_drutes $z $a $b 
      fi
    fi
    if [[  $l == "t"  ]]; then
      let z=$z+1
      if [[ $z -lt $nproc ]] ; then
        run_drutes $z $a $b $l  &
      else
        run_drutes $z $a $b $l
      fi
    fi
  done < pars.in
 

#at the end of drutes function file obj.val is created, if all files for each process exist then we have finished 
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

for i in `seq 1 $nproc` ; do
  val=`cat $i/objfnc.val | grep -v "#" `
  echo $i $val >> drutes.vals
done


  
