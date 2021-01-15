#!/bin/bash
#
# M. Leps  2006
#

if [ "$1" = "" ] ; then
    echo
    echo "Usage: doonperuns switch"
    echo
    exit 1
fi

list=`cat machines`

for i in $list ; do
  echo $i ;

# Make leps and osmina DIRECTORY
  if [ "$1" = "mkdir" ] ; then
    ssh $i " if [ ! -d /store/leps/ ] ; then mkdir /store/leps ; fi ; if [ ! -d /store/leps/osmina ] ; then mkdir /store/leps/osmina ; fi ;" ;
  fi

# COPY of DATA directory
  if [ "$1" = "cp" ] ; then
    ssh $i "rm -rf /store/leps/osmina/data ; cp -r /home/leps/cpp/oofem/osmina/data/ /store/leps/osmina/ " ;
  fi

# KILL ALL
  if [ "$1" = "kill" ] ; then
    ssh $i "killall much.sh davka.bat run_pasade.sh pasade run_oofem.sh oofem" ;
  fi
								
# Vytvori achiv
  if [ "$1" = "archiv" ] ; then
    ssh $i "cd /store/leps/osmina/data ; tar -czf pasade_$i.tgz uni-ns*.* hydro-ns*.* triax-ns*.* *.sav* *.log results* "
  fi

# PS to ALL
  if [ "$1" = "ps" ] ; then
    ssh $i "ps -u leps"
  fi

# RUN single JOB
  if [ "$1" = "run" ] ; then
    ssh $i "cd /store/leps/osmina/data ; ./much.sh 5 ./run_pasade.sh > pokus_$i.log 2>&1 & " &
  fi

# create directory for RESULTS
  if [ "$1" = "dir" ] ; then
    if [ ! -d "/home/leps/cpp/oofem/osmina/data/$i" ]; then
      echo "Creating directory /home/leps/cpp/oofem/osmina/data/$i " ;
      mkdir /home/leps/cpp/oofem/osmina/data/$i ;
    fi
  fi

# COPY RESULTS
  if [ "$1" = "store" ] ; then
    ssh $i "cd /store/leps/osmina/data ; cp results* *.log *.out *.tgz ~/cpp/oofem/osmina/data/$i/."
  fi

  if [ "$1" = "ls" ] ; then
    ssh $i "cd /store/leps/osmina/data ; ls -al *.out* "
  fi

  if [ "$1" = "ll" ] ; then
    ssh $i "cd /store/leps/osmina/data ; ls -al "
  fi

  if [ "$1" = "tgz" ] ; then
      ssh $i "cd /store/leps/osmina/data ; ls -al *.tgz"
  fi

  if [ "$1" = "log" ] ; then
        ssh $i "cd /store/leps/osmina/data ; cat *.log"
  fi

  if [ "$1" = "do" ] ; then
        ssh $i "$2 $3 $4 $5 $6"
  fi

done
	
exit 0

