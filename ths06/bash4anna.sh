#!/bin/bash

function run_drutes {

  echo "running id:" $2 "at node:" $1
 
  rm -rf $1

  mkdir $1

  cp -a drutemp/* $1

  cd $1
  
  a=$4
  n=$5
  Ks=$6
  ths="0.6"
  Ss="0.0"

  m=`echo "scale=12; 1.0-1.0/$n" | bc`
   
  sed -e 's/!A/'$a'/g' -e 's/!N/'$n'/g' -e 's/!M/'$m'/g' -e 's/!K/'$Ks'/g' -e 's/!T/'$ths'/g' -e 's/!S/'$Ss'/g' drutes.conf/water.conf/matrix_template.conf > drutes.conf/water.conf/matrix.conf
  
  bin/drutes -o optim > /dev/null
  
  
  cat out/obspt_RE_matrix-1.out | grep -v "#" > ../files/$3/$2
  

  touch objfnc.val

  cd ..

 
    
}


function dataprep {
      let z=0
      while read l idd a b c
        do
          if [[  $l == "p"  ]]; then
            let z=z+1
            if [[ $z -lt $nproc ]] ; then
              run_drutes $z $idd  $1 $a $b $c   &
            else
              run_drutes $z $idd $1 $a $b $c
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
 } 

function go {
 
    let nproc=25

    let count=0

    let id=0



    rm -f pars.in
    rm -rf files/$1
    mkdir files/$1
    while read  a b c 
      do
        let count=count+1
        let id=id+1
        echo "p"  $id $a $b $c  >> pars.in

        if [[  $count == $nproc  ]]; then
          dataprep $1
          echo "tam"
          let count=0 
          rm -rf pars.in
          echo $count $id
        fi
      done < inputs/$1
}
  

files[1]="X12Dgrid.txt"
files[2]="X2009.txt"


for item in ${files[*]}
do
    echo "processing file:" $item
    go $item
done
