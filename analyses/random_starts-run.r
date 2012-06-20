source("../src/r/multidimevol-lib.r")

# Require the argument --run=<runname> on the command line
the.args <- commandArgs()
run.name <- sub("--run=", "", the.args[grep("--run=", the.args)])
stopifnot(length(run.name) == 1)

# run.name for the submitted analyses:
# run.name <- "2011_10_20"

source(project.path("analyses/random_starts-setup.r"))

desired.runs <- 200
models <- list()
successes <- 0
runs <- lapply(1:nrow(choices.mx), function(i) {
   models[[i]] <- alter.model.list(m.start, as.list(choices.mx[i,]))
   if (abs.fitness(models[[i]]) <= 0) return(NULL)
   successes <<- successes + 1
   if (successes > desired.runs) return(NULL)
   cat(".")
   evolve(200, models[[i]], tolerance=1e-10, max.step=0.10)
})
cat("\n")

# only consider runs that start with fitness > 1. This is arbitrary, but prevents super big adaptive runs,
#  and if one accepts these as absolute fitnesses, then perhaps relative fitnesses can be interpretable
#  as selection coefficients.
runs <- runs[!sapply(runs, is.null)]
runs <- runs[sapply(runs, function(x) x[1,"fitness"]) > 0]

stopifnot(length(runs) >= desired.runs)
runs <- runs[1:desired.runs]

runs.collated <- collate.par.vs.time(runs)

ranges <- lapply(c("fitness", names(choices.mx)), function(param) range(unlist(lapply(runs, function(x) range(x[,param])))))
names(ranges) <- c("fitness", names(choices.mx))

save.image(file=project.path(paste("analyses/data/random_starts-", run.name, ".rimage", sep="")))

# END


