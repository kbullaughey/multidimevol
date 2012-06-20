#!/bin/bash
#$ -cwd
#$ -l h_vmem=2.5g
#$ -o grid/$JOB_ID.out
#$ -e grid/$JOB_ID.err

if [ "x" == "x$1" ] ; then
   echo "Must specify R script" >&2
   exit 1
fi
script=$1
shift

if [ ! -f $script ] ; then
   echo "File $script doesn't exist" >&2
   exit 2
fi

cmd="R --vanilla --args $* < $script"
echo $cmd
echo $cmd | bash

# END
