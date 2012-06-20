#!/bin/bash

if [ "x" == "x$1" ] ; then
  echo "Must specify run label" >&2
  exit 1
fi
run=$1

combos="data/dimension_analysis-combos-$run.txt"
if [ ! -f $combos ] ; then
  echo "Couldn't find combos file: $combos" >&2
  exit 1
fi

for i in `cat $combos`; do 
  rimage="data/dimension_analysis/$run/dimsim-$i.rimage"
  if [ ! -f $rimage ] ; then
    echo "Couldn't find rimage file: $rimage" >&2
    exit 1
  fi
  R --vanilla --args --rimage=$rimage < dimension_analysis-post_process.r
done
