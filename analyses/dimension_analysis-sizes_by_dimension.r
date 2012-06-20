# This script locates the reversals and computes the reversal delta and 
#  fitness increase for each, which are later used for a figure
library(grid)
source("../src/r/multidimevol-lib.r")

# Require the argument --run=<runname> on the command line
the.args <- commandArgs()
run.name <- sub("--run=", "", the.args[grep("--run=", the.args)])
stopifnot(length(run.name) == 1)

# base is the directory in which to find the *.rimage files
base <- project.path(paste("analyses/data/dimension_analysis/", run.name, sep=""))

# get the list of files
cmd <- paste("ls ", base, "/*biggest_reversals.rimage", sep="")
pip <- pipe(cmd, "r")
files <- scan(what="", file=pip, sep="\n")
close(pip)

# extract the parameter strings from the file names
param.strings <- sub("^.*/dimsim-(.*)-biggest_reversals.rimage$", "\\1", files)

# read in the reversals from the files
reversals <- lapply(files, function(f) {
   load(f)
   biggest.reversals
})
names(reversals) <- param.strings

at.least.one.reversal <- sapply(reversals, length)

param.list <- strsplit(param.strings, split="-")

bar.x <- 0.07
num.runs <- 31
runs.use <- 2
evolvable.traits <- c("Y.a", "Y.b", "Z.a", "Z.b", "k")
num.traits <- length(evolvable.traits)
param.expr <- c(expression(alpha[y]), expression(beta[y]), expression(alpha[z]), expression(beta[z]), "k")
small.plot.height <- 0.80
small.plot.width <- 0.90
horizontals <- c(0, 0.17, 0.75, 0.95)
heights <- (c(horizontals[-1], 1)-horizontals)[-3]
horizontals <- cumsum(c(0, heights/sum(heights)))[-4]
sep.width <- 0.3
set.sizes <- c(5,10,10,5,1)

# get binary grid indicating parameter presence according to our bitmask
combos <- t(expand.grid(rep(list(0:1), length(evolvable.traits))))
combos <- combos[,apply(combos, 2, sum)>0]

# Here I re-order the runs to get the pattern I want. 
combos.o <- order(apply(combos, 2, sum), 1-combos[1,], 1-combos[2,], 1-combos[3,], 1-combos[4,])
combos.ordered <- combos[,combos.o]

param.binary.strings <- sapply(param.list, function(x) paste((evolvable.traits %in% x)+0, collapse=""))
combos.binary.strings <- apply(combos.ordered, 2, paste, collapse="")

# reorder things
param.list.ordered <- param.list[match(combos.binary.strings, param.binary.strings)]
reversals.ordered <- reversals[match(combos.binary.strings, param.binary.strings)]
at.least.one.rev.ordered <- at.least.one.reversal[match(combos.binary.strings, param.binary.strings)]

# used for normalizing reversals
params <- rownames(reversals.ordered[[6]][[1]])
load(project.path(paste("analyses/data/mean_param_improvement-", run.name, ".rimage", sep="")))

run.dimensionality <- sapply(strsplit(sub("^.*dimsim-", "", names(reversals.ordered)), "-"), length)
reversal.sizes <- lapply(reversals.ordered, function(set) {
  if (length(set) == 0) return(numeric())
  sizes <- lapply(set, function(x) {
    unlist(x[,"reversal.delta"])/mean.param.improvement[1:5]
  })
  as.numeric(na.omit(unlist(sizes)))
})
reversal.sizes.by.dimension <- lapply(split(reversal.sizes, run.dimensionality), function(x) as.numeric(unlist(x)))

pdf(project.path(paste("analyses/plots/reversal_sizes_by_dimension-", run.name, ".pdf", sep="")), height=5, width=6)
par(mar=c(4,4,1,1), mgp=c(2.4,0.8,0), cex.axis=0.75)
boxplot(reversal.sizes.by.dimension, notch=TRUE, outcex=0.7, col="gray80",
  outcol="gray50",
  xlab="Number of evolvable parameters", ylab=expression(paste(plain("Normalized reversal size "), (Delta))))
text(x=1, y=c(0.09, 0.04), labels=c("no", "reversals"))
dev.off()


# END
