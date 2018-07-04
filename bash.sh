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

  
#   a=`echo ${ka} | sed -e 's/[eE]+*/\\*10\\^/'`
#   kd=`echo ${kd} | sed -e 's/[eE]+*/\\*10\\^/'`   
# 
#   kaf=`echo ${kaf} | sed -e 's/[eE]+*/\\*10\\^/'`
#   Ks=`echo ${Ks} | sed -e 's/[eE]+*/\\*10\\^/'`



  
#  
#   ka=`echo "e($ka)" | bc -l`
# #  kd=`echo "e($kd)" | bc -l`
#  
#   kaf=`echo "e($kaf)" | bc -l`
#   kdf=`echo "e($kdf)" | bc -l`

m=`echo "1-1/$n" | bc -l`

 
   
#   substitution of parameters into input files for drutes 
  sed -e 's/!a/'$a'/g' -e 's/!n/'$n'/g'  -e 's/!m/'$m'/g' -e 's/!ths/'$ths'/g'  -e 's/!Ks/'$Ks'/g'  drutes.conf/water.conf/matrix.conf.temp > drutes.conf/water.conf/matrix.conf

# sed -e 's/!A/'$a'/g' -e 's/!N/'$n'/g' -e 's/!M/'$m'/g' -e 's/!K/'$Ks'/g' -e 's/!T/'$ths'/g' -e 's/!S/'$Ss'/g' drutes.conf/water.conf/matrix.conf.temp > drutes.conf/water.conf/matrix.conf

  
  #execute drutes, send the output into "black hole"
  bin/drutes -o optim --tmax 17 s > /dev/null
#   
#   let val=0
#   
  cat out/objfnc.val | grep -v '#' > objfnc.val
#   
#   while read n
#     do   
#       n=`echo ${n} | sed -e 's/[eE]+*/\\*10\\^/'`
#       val=`echo "$val+$n" | bc -l`
#     done < ll
#     
#   
#   echo $val > objfnc.val
#   echo " " >> objfnc.val
#   
#   gnuplot < ../plot.gnuplot
  
  echo $a $n $ths $Ks  $val >> ../totvals
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
while read l a b c d
  do
    if [[  $l == "p"  ]]; then
      let z=z+1
      if [[ $z -lt $nproc ]] ; then
        run_drutes $z $a $b $c $d  &
      else
        run_drutes $z $a $b  $c $d
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


  
