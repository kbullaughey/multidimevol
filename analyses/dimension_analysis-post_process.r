# This script locates the reversals and computes the reversal delta and 
#  fitness increase for each, which are later used for a figure

source("../src/r/multidimevol-lib.r")
the.args <- commandArgs()
rimage <- sub("--rimage=", "", the.args[grep("--rimage=", the.args)])
load(rimage)

param.str <- sub("\\.rimage", "", sub(".*dimsim-", "", rimage))
evolving.traits <- strsplit(param.str, "-")[[1]]

ranges <- lapply(c("fitness", names(choices.mx)), function(param) range(unlist(lapply(runs, function(x) range(x[,param])))))
names(ranges) <- c("fitness", names(choices.mx))

process.tick <- function(tick, run, param) {
   p <- run[[param]] > tick
   crosses <- which(p[-1] != p[-length(p)])
   if (length(crosses) < 2) return(c(reversal.delta=NA, fitness.delta=NA, fitness2.delta=NA, peak=NA))
   this.res <- lapply(2:length(crosses), function(i) {
      if (crosses[i]-1 == crosses[i-1])
         return(NULL)
      cross.indicies <- crosses[i-1]:crosses[i]
      reversal.delta <- (function(x) x[2]-x[1])(range(run[[param]][cross.indicies]))
      fitness.delta <- run$fitness[crosses[i]] - run$fitness[crosses[i-1]]
      wmax <- which.max(run[[param]][cross.indicies])
      wmin <- which.min(run[[param]][cross.indicies])
      peak.i <- cross.indicies[c(wmax,wmin)[which.max((mean(run[[param]][crosses[(i-1):i]])-run[[param]][cross.indicies][c(wmax,wmin)])^2)]]
#      peak <- c(wmax,wmin)[which.max((mean(run[[param]][crosses[(i-1):i]])-run[[param]][c(wmax,wmin)])^2)]
#      peak.i <- (crosses[i-1]:crosses[i])[peak]
      second.half.fitness.delta <- run[["fitness"]][crosses[i]] / run[["fitness"]][peak.i] - 1
      return(c(reversal.delta=reversal.delta, fitness.delta=fitness.delta, 
         fitness2.delta=second.half.fitness.delta, peak=peak.i))
   })
   this.res <- this.res[!sapply(this.res, is.null)]
   if (length(this.res) > 0) {
      return(this.res[[which.max(sapply(this.res, function(x) x[1]))]])
   } else {
      return(c(reversal.delta=NA, fitness.delta=NA, fitness2.delta=NA, peak=NA))
   }
}

params <- NULL
biggest.reversals <- lapply(runs, function(run) {
   run <- as.data.frame(run)
   params <<- names(run)[-(1:2)]
   t(sapply(params, function(param) {
      ticks <- (function(x) seq(x[1], x[2], length.out=40))(range(run[[param]]))
      res <- as.data.frame(t(sapply(ticks, process.tick, run, param)))
      if (sum(!is.na(res$rev)) == 0) return(c(reversal.delta=NA, fitness.delta=NA, fitness2.delta=NA, peak=NA))
      return(res[which.max(res[,1]),])
   }))
})
biggest.reversals <- biggest.reversals[sapply(biggest.reversals, function(x) sum(!is.na(x)))>0]
out.file <- sub("\\.rimage", "-biggest_reversals.rimage", rimage)
save(biggest.reversals, file=out.file)

# END




# END
