source("../src/r/multidimevol-lib.r")

# Require the argument --run=<runname> on the command line
the.args <- commandArgs()
run.name <- sub("--run=", "", the.args[grep("--run=", the.args)])
stopifnot(length(run.name) == 1)

# This script sets up the 31 combinations of parameters I want to test.
# These sets are all the possible non-empty subsets of the evolvable traits:
evolvable.traits <- c("Y.a", "Y.b", "Z.a", "Z.b", "k")

combos <- expand.grid(rep(list(0:1), length(evolvable.traits)))
combos <- combos[apply(combos, 1, sum)>0,]

combos.et <- apply(combos, 1, function(x) evolvable.traits[x==1])
runs <- paste("qsub r_submit.sh dimension_analysis-run.r --run=", run.name, 
  " --evolvable=", sapply(combos.et, paste, collapse=","), sep="")

shell.script <- project.path(paste("analyses/dimension_analysis-launch-", run.name, ".sh", sep=""))
cat("#!/bin/bash\n", file=shell.script)
cat(file=shell.script, runs, sep="\n", append=TRUE)

# END
