#!/bin/bash

function run_drutes {

  rm -rf $1

  mkdir $1

  cp -a drutemp/* $1

  cd $1
  
  
  
  ka=`echo "e($2)" | bc -l`
  kd=$3
 # kd=`echo "e($3) | bc -l"`
  



  sed -e 's/!ka/'$ka'/g' -e 's/!cmax/'$kd'/g'  drutes.conf/ADE/sorption_temp.conf > drutes.conf/ADE/sorption.conf
  
  #execute drutes, send the output into "black hole"
  bin/drutes > /dev/null
  
  let val=0
  
  cat out/objfnc.val | grep -v '#' > ll
  
  while read n
    do   
      n=${n/[eE]+*/*10^}
      val=`echo "$val+$n" | bc -l`
    done < ll
    
  
  echo $val > objfnc.val
  echo " " >> objfnc.val

  cd ..
    
}


rm -f drutes.vals

#count the number of processes       
let nproc=0
while read l a b  ; do
    if [[  $l == "p"  ]]; then
      let nproc=nproc+1
    fi
  done < pars.in
  
#execute drutes function in parallel
let z=0
while read l a b
  do
    if [[  $l == "p"  ]]; then
      let z=z+1
      if [[ $z -lt $nproc ]] ; then
        run_drutes $z $a $b  &
      else
        run_drutes $z $a $b  
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


  
