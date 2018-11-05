#!/bin/bash
particle1=$1
particle2=$2

# for trunclevel in {0..95..5}
for ((trunclevel=0; trunclevel<100; trunclevel+=5))
do
  for ((trunclevel2=0; trunclevel2<100; trunclevel2+=5))
  do
    if (($trunclevel + $trunclevel2 > 90))
    then
      continue
    fi
    printf "**************************************  $trunclevel , $trunclevel2  **************************\n"
    root -l -b -q p2pisep.C+"($trunclevel,$trunclevel2,\"$particle1\",\"$particle2\")"
  done
done
