library(grid)
library(xtable)
source("../src/r/multidimevol-lib.r")

# Require the argument --run=<runname> on the command line
the.args <- commandArgs()
run.name <- sub("--run=", "", the.args[grep("--run=", the.args)])
stopifnot(length(run.name) == 1)

# run.name for the submitted analyses:
# run.name <- "2011_10_20"

shuffle <- function(x, len=length(x)) x[sample(1:length(x), len, replace=FALSE)]

reduce <- function(x, len=100) x[floor(seq(1, length(x), length.out=len))]
thin <- function(x, step=0.01) {
  if (length(x) == 0) returnr(c())
  use <- rep(FALSE, length(x))
  use[1] <- TRUE
  if (length(x) == 1) return(use)
  prev <- x[1]
  for (i in 2:length(x)) {
    if (x[i] >= prev+step) {
      use[i] <- TRUE
      prev <- x[i]
    }
  }
  return(use)
}

load(file=project.path(paste("analyses/data/random_starts-", run.name, ".rimage", sep="")))
num.runs <- length(runs)

params <- colnames(runs[[1]])[-(1:2)]
time.range <- range(sapply(runs, function(x) range(x[,"time"])))
param.ranges <- sapply(params, function(p) expand.range2(range(sapply(runs, function(x) range(x[,p])))))

min.fit <- min(sapply(runs, function(x) min(x[,"fitness"])))
max.fit <- max(sapply(runs, function(x) max(x[,"fitness"])))
fit.step.size <- (max.fit - min.fit)/1000

runs.reduced <- lapply(runs, function(r) {
   cat(".")
   r[thin(r[,"fitness"], step=fit.step.size),]
})
cat("\n")

# I want something sensible to normalize the reversal delta by so that the reversals affecting 
#   different parameters are comparable, and I want it to be the same for all runs. 
# So I use as the denominator, the average change in the parameter over the course of the run
param.norm.factor <- sapply(params, function(p) 
  mean(sapply(runs, function(r) (function(x) abs(x[1]-x[length(x)]))(r[,p]))))

# determine the smaller of either the left or right trip to the peak
peak.reversal.delta <- function(before, peak, after, x) {
  res <- list(left=list(), right=list(), switch=peak)
  if (peak > 1)
    m <- max(x[peak] - x[before:peak])
    w <- which.max(x[peak] - x[before:peak])
    res$left <- list(max=m, which=w)
  if (peak < length(x))
    m <- max(x[peak] - x[peak:after])
    w <- which.max(x[peak] - x[peak:after]) + peak - 1
    res$right <- list(max=m, which=w)
  return(res)
}

# determine the smaller of either the left or right trip to the trough
trough.reversal.delta <- function(before, trough, after, x) {
  res <- list(left=list(), right=list(), switch=trough)
  if (trough > 1)
    m <- max(x[before:trough] - x[trough])
    w <- which.max(x[before:trough] - x[trough])
    res$left <- list(max=m, which=w)
  if (trough < length(x))
    m <- max(x[trough:after] - x[trough])
    w <- which.max(x[trough:after] - x[trough]) + trough - 1
    res$right <- list(max=m, which=w)
  return(res)
}

process.parameter <- function(run, param) {
  x <- run[[param]]
  left <- x[-c(length(x)-1, length(x))]
  middle <- x[-c(1,length(x))]
  right <- x[-(1:2)]
  peaks <- which(left < middle & middle > right)
  troughs <- which(left > middle & middle < right)
  boundaries <- sort(c(1,nrow(run),peaks,troughs))

  if (length(boundaries) <= 2) 
    return(NULL)

  peaks.before <- boundaries[match(peaks, boundaries)-1]
  peaks.after <- boundaries[match(peaks, boundaries)+1]
  troughs.before <- boundaries[match(troughs, boundaries)-1]
  troughs.after <- boundaries[match(troughs, boundaries)+1]
  
  reversal.delta.candidates <- list()
  if (length(peaks) > 0) {
    for (i in 1:length(peaks)) {
      reversal.delta.candidates[[length(reversal.delta.candidates)+1]] <- 
        peak.reversal.delta(peaks.before[i], peaks[i], peaks.after[i], x)
    }
  }
  if (length(troughs) > 0) {
    for (i in 1:length(troughs)) {
      reversal.delta.candidates[[length(reversal.delta.candidates)+1]] <- 
        trough.reversal.delta(troughs.before[i], troughs[i], troughs.after[i], x)
    }
  }
  
  # add on the fitness deltas
  reversal.delta.candidates <- lapply(reversal.delta.candidates, function(candidate) {
    if (length(candidate$left) > 0)
      candidate$left$fit.delta <- run$fitness[candidate$switch] - run$fitness[candidate$left$which]
      candidate$left$fit.start <- run$fitness[candidate$left$which]
    if (length(candidate$right) > 0)
      candidate$right$fit.delta <- run$fitness[candidate$right$which] - run$fitness[candidate$switch]
      candidate$right$fit.start <- run$fitness[candidate$switch]
    return(candidate)
  })

  return(reversal.delta.candidates)
}

filter.candidates <- function(reversal.delta.candidates) {
  rdc.i <- which.max(sapply(reversal.delta.candidates, function(candidate) min(candidate$left$max, candidate$right$max)))
  if (length(rdc.i) == 0) return(NULL)
  return(reversal.delta.candidates[[rdc.i]])
}

reformat.reversal <- function(candidate) {
  # use as the reversal delta, the lesser of the left and right distance
  # use as the fitness delta, the right fitness increase
  if (is.null(candidate)) 
    return(list(reversal.delta=NA, fitness.delta=NA, peak=NA))
  list(
    reversal.delta=min(candidate$left$max, candidate$right$max), 
    fitness.delta=candidate$right$fit.delta, 
    peak=candidate$switch
  )
}

mean.param.improvement <- sapply(c(params, "fitness"), function(param) {
  improvement <- mean(as.numeric(sapply(runs, function(x) max(x[,param])-min(x[,param]))))
})

biggest.reversals.stage1 <- lapply(runs, function(run) {
   run <- as.data.frame(run)
   lapply(params, function(param) process.parameter(run, param))
})

# examine the distribution of reversal sizes relative to the mean for each parameter and together
reversal.sizes <- lapply(biggest.reversals.stage1, function(x) {
  lapply(x, function(y) {
    sapply(y, function(z) min(z$left$max, z$right$max))
  })
})
# rearrang the list so that it's grouped by parameter
reversal.sizes.normalized <- lapply(1:length(params), function(i) {
  lapply(reversal.sizes, function(x) {
    if (length(x[[i]]) > 0) {
      res <- x[[i]]/mean.param.improvement[params[i]]
      names(res) <- NULL
      res
    } else  {
      numeric(0)
    }
  })
})
reversal.sizes.normalized.unlisted <- lapply(reversal.sizes.normalized, unlist) 

normalized.threshold <- 0.05

biggest.reversals.stage2 <- lapply(1:length(biggest.reversals.stage1), function(rep) {
  lapply(1:length(params), function(param.i) {
    param <- params[param.i]
    y <- biggest.reversals.stage1[[rep]][[param.i]]
    if (is.null(y)) return(list())
    sizes <- reversal.sizes.normalized[[param.i]][[rep]]
    # remove reversals that are too small
    stopifnot(length(y) == length(sizes))
    y[sizes > normalized.threshold]
  })
})

biggest.reversals.stage3 <- lapply(biggest.reversals.stage2, function(x) {
  # loop over the 5 parameters
  # select the best reversal, if there are more than one
  lapply(x, filter.candidates)
})

biggest.reversals.stage4 <- sapply(biggest.reversals.stage3, function(x) {
  # loop over the 5 parameters
  lapply(x, function(y) {
    reformat.reversal(y)
  })
})
rownames(biggest.reversals.stage4) <- params

# put together a matrix of the ones we want to display
biggest.reversals.mx <- sapply(params, function(p) {
  sapply(biggest.reversals.stage4[p,], function(x) x$reversal.delta)
})
fitness.delta.mx <- sapply(params, function(p) {
  sapply(biggest.reversals.stage4[p,], function(x) x$fitness.delta)
})

too.small <- t(t(biggest.reversals.mx) <= 0.05*mean.param.improvement[params])  |
  fitness.delta.mx <= 0.05*mean.param.improvement["fitness"]

biggest.reversals.mx[too.small] <- NA
fitness.delta.mx[too.small] <- NA

no.reversals.all <- which(apply(biggest.reversals.mx, 1, function(x) sum(!is.na(x))==0))

stopifnot(sum(params != c("Y.a", "Y.b", "Z.a", "Z.b", "k"))==0)
param.expr <- c(expression(alpha[y]), expression(beta[y]), 
   expression(alpha[z]), expression(beta[z]), expression(k))

set.seed(7)
# find the largest.bunch reversals to highlight
largest.bunch <- lapply(params, function(p) {
  x <- which(rank(biggest.reversals.mx[,p], na.last="keep") > sum(!is.na(biggest.reversals.mx[,p]))-10)
  if (length(x) > 5) {
    x <- shuffle(x, len=5)
  }
  x
})
names(largest.bunch) <- params

plot.which.col <- 1:5
plot.which.row <- rep(1, 5)
#param.axis.labels <- apply(param.ranges, 2, pretty, n=2)
param.axis.labels <- list(
  Y.a=c(0.4, 1.4, 2.4, 3.4),
  Y.b=c(0.2, 0.7, 1.2, 1.7),
  Z.a=c(0.6, 1.0, 1.4, 1.8),
  Z.b=c(25, 40, 55, 70),
  k=c(0, 0.2, 0.4, 0.6)
)
stopifnot(sum(names(param.axis.labels) != colnames(param.ranges)) == 0)
fit.axis.labs <- pretty(c(0, max.fit), n=3)
rev.counts <- apply(biggest.reversals.mx, 2, function(x) sum(x>0.05, na.rm=TRUE))
main.height <- 0.82
top.mar <- 0.1
left.mar <- 0.1

num.rev.reps <- 6
num.nonrev.reps <- 3
num.reps <- num.rev.reps + num.nonrev.reps
to.use <- as.numeric(unlist(largest.bunch))
to.use <- shuffle(to.use, len=num.rev.reps)
no.reversals.use <- intersect(no.reversals.all, which(sapply(runs.reduced, function(x) x[1,"fitness"]) < 4))
to.use <- c(to.use, head(no.reversals.use, num.nonrev.reps))
plot.cols <- c(1:num.rev.reps, (1:num.nonrev.reps)+num.rev.reps+1)
to.use[1] <- 5

scale <- 1.1
palette(c("black", rgb(0, 0, 1, 0.3), "coral3", "gray40", "gray70", "black"))
pdf(file=project.path(paste("analyses/plots/curves-grid-", run.name, "-v2.pdf", sep="")), 
  width=(0.2+num.reps)*scale, height=5.5*scale)
# full device
pushViewport(viewport())
  grid.text("evolvable traits", y=(2*top.mar+main.height)/2, x=0.02, rot=90, gp=gpar(cex=1.1))

  # top legend
  pushViewport(viewport(y=1, x=0.24, height=top.mar, width=0.25, just=c(0,1)))
    pushViewport(viewport(height=0.5))
      legend.y <- seq(1,0,len=3)
      grid.text(c("adaptive reversal", "trait does not show reversal", "other trajectories"), 
        y=legend.y, x=0.25, gp=gpar(cex=0.90), just="left")
      grid.segments(0.1, legend.y, 0.2, legend.y, gp=gpar(col=c(3, 6, 4), lwd=c(2,2,1)))
    popViewport()
  popViewport()
  
  # central main plot
  pushViewport(viewport(x=left.mar, y=1-top.mar, width=1-left.mar-0.02, height=main.height, just=c(0,1),
      layout=grid.layout(ncol=1, nrow=5+4, heights=rep(c(1,0.17), 5)[1:9])))
    grid.rect(gp=gpar(fill=5, col=5))
    for (i in 1:5) {
      p <- params[i]
      pushViewport(viewport(layout.pos.col=1, layout.pos.row=i*2-1))
        pushViewport(viewport(layout=grid.layout(widths=c(rep(1,num.rev.reps), 0.1, rep(1,num.nonrev.reps)), 
            ncol=num.reps+1, nrow=1)))
          grid.rect(gp=gpar(col=NA, fill="white"))
          pushViewport(viewport(layout.pos.row=1, layout.pos.col=num.rev.reps+1))
            grid.rect(gp=gpar(fill=5, col=5))
          popViewport()
          for (j in 1:num.reps) {
            w <- to.use[j]
            pushViewport(viewport(layout.pos.row=1, layout.pos.col=plot.cols[j], 
                xscale=range(fit.axis.labs), yscale=range(param.axis.labels[[i]])))
              grid.rect()
              grid.segments(unit(0, "npc"), unit(param.axis.labels[[i]], "native"), 
                unit(1, "npc"), unit(param.axis.labels[[i]], "native"), gp=gpar(col="gray80"))
              grid.segments(unit(fit.axis.labs, "native"), 
                unit(1, "npc"), unit(fit.axis.labs, "native"), unit(0, "npc"), gp=gpar(col="gray80"))
              if (j==1) {
                grid.yaxis(at=param.axis.labels[[i]], gp=gpar(cex=0.80))
                grid.text(param.expr[i], y=0.5, x=-0.54, gp=gpar(cex=1.3))
                if (i == 1) {
                  grid.xaxis(main=FALSE, at=fit.axis.labs, gp=gpar(cex=0.80))
                  grid.text("fitness", x=0.5, y=1.54, gp=gpar(cex=1.0))
                }
              }
              if (j <= num.rev.reps) {
                subset <- setdiff(1:nrow(biggest.reversals.mx), no.reversals.all)
              } else {
                subset <- no.reversals.all
              }
              if (sum(j == c(1,7))) {
                for (k in 1:15) {
                  grid.lines(x=reduce(runs.reduced[subset][[k]][,"fitness"]), y=reduce(runs.reduced[subset][[k]][,p]), 
                    default.units="native", gp=gpar(col=2))
                }
              }
              is.reversal <- !is.na(biggest.reversals.mx[w,p])
              grid.lines(x=runs.reduced[[w]][,"fitness"], y=runs.reduced[[w]][,p], 
                default.units="native", gp=gpar(col=c(6,3)[1+is.reversal], lwd=c(2,2)[1+is.reversal]))
            popViewport()
          }
          if (i==5) {
            pushViewport(viewport(layout.pos.row=1, layout.pos.col=c(1,num.rev.reps), xscale=c(0,num.rev.reps)))
              grid.text("at least one adaptive reversal", x=0.5, y=-0.4, gp=gpar(cex=1.0))
              grid.text("replicate:", x=-0.00, y=-0.15, gp=gpar(cex=0.9), just="right")
              grid.text(1:num.rev.reps, x=unit((1:num.rev.reps)-0.5, "native"), y=-0.15, gp=gpar(cex=0.9), just="right")
            popViewport()
            pushViewport(viewport(layout.pos.row=1, layout.pos.col=num.rev.reps+1+c(1,num.nonrev.reps), xscale=c(0,num.nonrev.reps)))
              grid.text("no adaptive reversals", x=0.5, y=-0.4, gp=gpar(cex=1.1))
              grid.text((1:num.nonrev.reps)+num.rev.reps, x=unit((1:num.nonrev.reps)-0.5, "native"), y=-0.15, 
                gp=gpar(cex=0.9), just="right")
            popViewport()
          }
        popViewport()
      popViewport()
    }
  popViewport()
popViewport()
dev.off()

save.image(project.path(paste("analyses/data/post-random_starts-curves-grid-", run.name, ".rimage", sep="")))
save(mean.param.improvement, file=project.path(paste("analyses/data/mean_param_improvement-", run.name, ".rimage", sep="")))

# generate a table similar to Table 2, but for the SMS approximation

# right reversal sizes
right.reversal.deltas <- lapply(1:length(params), function(i) {
  unlist(lapply(biggest.reversals.stage3, function(x) 
    as.numeric(x[[i]]$right$max / mean.param.improvement[params[i]])))
})
right.medians <- sapply(right.reversal.deltas, median)
right.rel.fitness.deltas <- lapply(1:length(params), function(i) {
  unlist(lapply(biggest.reversals.stage3, function(x)
    x[[i]]$right$fit.delta / x[[i]]$right$fit.start))
})
right.rel.fitness.deltats.medians <- sapply(right.rel.fitness.deltas, median)

# left reversal sizes
left.reversal.deltas <- lapply(1:length(params), function(i) {
  unlist(lapply(biggest.reversals.stage3, function(x) 
    as.numeric(x[[i]]$left$max / mean.param.improvement[params[i]])))
})
left.medians <- sapply(left.reversal.deltas, median)
left.rel.fitness.deltas <- lapply(1:length(params), function(i) {
  unlist(lapply(biggest.reversals.stage3, function(x)
    x[[i]]$left$fit.delta / x[[i]]$left$fit.start))
})
left.rel.fitness.deltats.medians <- sapply(left.rel.fitness.deltas, median)

rev.counts <- sapply(right.rel.fitness.deltas, length)

reversal.table <- data.frame(
  Trait=params, 
  Reversals=rev.counts, 
  for.med.s=left.rel.fitness.deltats.medians,
  for.med.delta=left.medians,
  back.med.s=right.rel.fitness.deltats.medians,
  back.med.delta=right.medians
)

reversal.xtable <- xtable(reversal.table)
print(reversal.xtable, file=project.path(paste("analyses/data/SMS_reversals.table-", run.name, ".tex", sep="")))

# END



