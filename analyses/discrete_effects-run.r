source("../src/r/multidimevol-lib.r")

the.args <- commandArgs()
run <- sub("--run=", "", the.args[grep("--run=", the.args)])
set <- sub("--set=", "", the.args[grep("--set=", the.args)])
regconfig <- sub("--regconfig=", "", the.args[grep("--regconfig=", the.args)])
simconfig <- sub("--simconfig=", "", the.args[grep("--simconfig=", the.args)])

stopifnot(length(run)==1)
stopifnot(length(regconfig)==1)
stopifnot(length(simconfig)==1)

# run a simulation 
reg.mod.start <- load.config(regconfig, defaults=default.model)
sim.mod.start <- load.config(simconfig, defaults=default.simulation)
run.res <- popsim(mod.reg=reg.mod.start, mod.sim=sim.mod.start)

out.dir <- project.path(paste("analyses/data/discrete_effects-simulations/", set, "/", run, sep=""))
if (!file.exists(out.dir)) {
   dir.create(out.dir, recursive=TRUE)
}

out <- t(sapply(run.res$fixations, unlist))
write.config.file(paste(out.dir, "/", set, "-", run, "-before.rconf", sep=""), reg.mod.start)
write.config.file(paste(out.dir, "/", set, "-", run, ".sconf", sep=""), sim.mod.start)
write.table(out, file=paste(out.dir, "/", set, "-", run, ".out", sep=""), quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")
save(run.res, reg.mod.start, sim.mod.start, file=paste(out.dir, "/", set, "-", run, ".rimage", sep=""))
write.config.file(paste(out.dir, "/", set, "-", run, "-after.rconf", sep=""), run.res$mod.final)

# END



