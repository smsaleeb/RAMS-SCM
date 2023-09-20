#!/bin/bash

rt=$1
otherversion1=$2
otherversion2=$3

if [ ! -n "$rt" ]; then
 echo "Need command line argument for RUNTYPE"
 echo "1 = Check scm.in   vs  scm.out"
 echo "2 = Check Argument-2 directory  vs  Argument-3 directory"
 exit
fi

if [ $rt -eq 1 ]; then
 for f in scm.in/*; do
  file=`basename "$f"`
  echo "DIFF: $file "
  diff $f scm.out/$file
 done
fi

if [ $rt -eq 2 ]; then
 if [ -n "$otherversion1" -a -n "$otherversion2" ]; then
  for f in $otherversion1/*; do
   file=`basename "$f"`
   echo "DIFF: $file "
   if [ -f $otherversion2/$file ]; then diff $f $otherversion2/$file; fi
  done
 fi
fi
