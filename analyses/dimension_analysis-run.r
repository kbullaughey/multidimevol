source("../src/r/multidimevol-lib.r")

the.args <- commandArgs()
evolvable <- sub("--evolvable=", "", the.args[grep("--evolvable=", the.args)])
run.name <- sub("--run=", "", the.args[grep("--run=", the.args)])

all.params <- c("Y.a", "Y.b", "Z.a", "Z.b", "k")

# make sure we have a new place to put the runs
base.dir <- project.path(paste("analyses/data/dimension_analysis/", run.name, sep=""))
if (!file.exists(base.dir)) {
   if (!dir.create(base.dir, recursive=TRUE)) {
      cat(file=stderr(), "Exiting...failed to create directory\n")
      quit()
   }
}

# check and parse parameters
stopifnot(length(evolvable) == 1)
evolvable <- strsplit(evolvable, split=",")[[1]]

# compute the bitmask for this selection of parameters
ep.bitmask <- sum(2^(match(evolvable, all.params)-1))

# this provides:
#  m.start
#  choices.mx
source("random_starts-setup.r")

models <- list()
runs <- lapply(1:1000, function(i) {
   models[[i]] <<- alter.model.list(m.start, as.list(choices.mx[i,][evolvable]))
   full <- evolve(200, models[[i]], tolerance=1e-10, max.step=0.10, evolvable.param.bitmask=ep.bitmask)
   # Thin the data a bit so file sizes aren't so large
   full[seq(1, nrow(full), by=8),]
})

# only consider runs that start with fitness > 0. This is arbitrary, but prevents super big adaptive runs,
#  and if one accepts these as absolute fitnesses, then perhaps relative fitnesses can be interpretable
#  as selection coefficients.
runs <- runs[sapply(runs, function(x) x[1,"fitness"]) > 0]

if (length(runs) < 200) {
   cat(file=stderr(), "Failed to start with 200 random positive fitnesses\n")
   quit()
}
runs <- runs[1:200]

runs.collated <- collate.par.vs.time(runs)

save.image(file=paste(base.dir, "/dimsim-", paste(evolvable, collapse="-"), ".rimage", sep=""))


# END



