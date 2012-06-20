#!/bin/bash
#$ -l h_vmem=1g
#$ -o grid/closed.$JOB_ID.out
#$ -e grid/closed.$JOB_ID.err
#$ -cwd

rscript="discrete_effects-run.r"

if [ "x$1" == "x" ] ; then
   echo "Must give run name" >&2
   exit 1
fi
run=$1

if [ "x$2" != "x" ] ; then
   regconfig=$2
else
   regconfig="configs/${run}.rconf"
fi
if [ ! -f $regconfig ] ; then
   echo "Failed to find regulatory model file: $regconfig" >&2
   exit 1
fi

if [ "x$3" != "x" ] ; then
   simconfig=$3
else
   simconfig="configs/${run}.sconf"
fi
if [ ! -f $simconfig ] ; then
   echo "Failed to find simulation model file: $simconfig" >&2
   exit 1
fi

R --vanilla --args --run=$run-$JOB_ID --prefix=runs --regconfig=$regconfig \
  --simconfig=$simconfig < $rscript

R --vanilla --args --rimage=runs/$run-$JOB_ID/$run-$JOB_ID.rimage \
  --pdf=runs/$run-$JOB_ID/$run-$JOB_ID-selection.pdf < selection-vertical.r

# END
