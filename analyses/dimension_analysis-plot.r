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

max.rev <- max(unlist(lapply(reversals, function(x) {
   lapply(x, function(y) {
      y[,"reversal.delta"]
   })
})), na.rm=TRUE)
max.fi <- max(unlist(lapply(reversals, function(x) {
   lapply(x, function(y) {
      y[,"fitness2.delta"]
   })
})), na.rm=TRUE)

run.rev.count <- numeric(num.runs)
palette(c("#CD5B45", "#BF801B", "#41A521", "#1874CD", "#9A69AF"))
pdf(file=project.path(paste("analyses/plots/dimension_analysis-", run.name, ".pdf", sep="")), height=3.00, width=8)
# main viewport
pushViewport(viewport())
  # middle viewport
  pushViewport(viewport(x=0.5, y=mean(horizontals[2:3]), height=(horizontals[3] - horizontals[2])*0.95, width=1))
    grid.segments(bar.x, 0, bar.x, 1)
    # left labels viewport
    pushViewport(viewport(x=bar.x/2, y=0.5, width=bar.x, height=1, yscale=c(0,num.traits)))
      grid.text(param.expr, x=0.7, y=unit((num.traits:1)-0.5, "native"), gp=gpar(cex=1.2))
      grid.text("evolvable trait", x=0.2, y=0.5, rot=90, gp=gpar(cex=0.9))
    popViewport()
    # right main viewport
    col.widths <- c(rep(1, 5), sep.width, rep(1, 10), sep.width, rep(1,10), sep.width, rep(1,5), sep.width, 1)
    pushViewport(viewport(x=(1+bar.x)/2, y=0.5, height=1, width=(1-bar.x)*0.99, 
         layout=grid.layout(nrow=num.traits, ncol=num.runs+num.traits-1, widths=col.widths)))
      grid.rect(gp=gpar(col="gray90", fill="gray90"))
#      grid.segments(0, (1:4)/5, 1, (1:4)/5, gp=gpar(col="gray70"))
      # grid of small viewports
      prev.num <- 1
      for (run in 1:num.runs) {
         this.num.traits <- sum(combos.ordered[,run])
         if (prev.num < this.num.traits) {
            pushViewport(viewport(layout.pos.row=c(1,num.traits), layout.pos.col=run+this.num.traits-2))
               grid.segments(0.5, 0, 0.5, 1, gp=gpar(col="gray50"))
            popViewport()
         }
         prev.num <- this.num.traits
         for (trait.i in 1:num.traits) {
            if (combos.ordered[trait.i,run]) {
               pushViewport(viewport(layout.pos.row=trait.i, layout.pos.col=run+this.num.traits-1))
                  reversal.delta <- sapply(reversals.ordered[[run]], function(x) x[trait.i,"reversal.delta"][[1]])
                  fit.improve <- sapply(reversals.ordered[[run]], function(x) x[trait.i,"fitness2.delta"][[1]])
                  grid.rect(gp=gpar(col="gray90", fill="white"))
                  pushViewport(viewport(height=small.plot.height, width=small.plot.width, xscale=c(0,max.rev), yscale=c(0,max.fi)))
                     rev.count <- round(sum(!is.na(reversal.delta))/2)
                     if (rev.count > 0) {
                        grid.text(rev.count, x=0.5, y=0.5, gp=gpar(cex=0.65, col="#CD5B45"))
                        run.rev.count[run] <- run.rev.count[run] + rev.count
                     } else {
                        grid.text(rev.count, x=0.5, y=0.5, gp=gpar(cex=0.65, col="#1874CD"))
                     }
                  popViewport()
               popViewport()
            }
         }
      }
    popViewport()
  popViewport()
  # push the bottom viewport
  pushViewport(viewport(x=(bar.x+1)/2, y=mean(horizontals[1:2]), height=horizontals[2]-horizontals[1], width=(1-bar.x)*0.99, 
     xscale=c(0,num.runs+num.traits-1), yscale=0:1))
    grid.segments(c(0,6,17,28,34),0.9,c(5,16,27,33,35),0.9, default.units="native", gp=gpar(col=1:5, lwd=2))
    grid.text(paste(1:5, c(rep(" trait",4), ""), c("", rep("s", 3), ""), sep=""), 
      x=unit(c(2.5, (6+16)/2, (17+27)/2, (28+33)/2, 34.5), "native"), y=0.7, gp=gpar(cex=1.0, col=1:5))
    grid.text("dimensionality", x=0.5, y=0.25, gp=gpar(cex=1.2))
  popViewport()
  # push the top viewport
  pushViewport(viewport(x=(bar.x+1)/2, y=mean(horizontals[3:4]), height=horizontals[4]-horizontals[3], width=(1-bar.x))) 
    mean.rev <- run.rev.count/rep(1:5, set.sizes)
    at.least.one.rev.ordered.percent <- at.least.one.rev.ordered/200*100
    pushViewport(viewport(x=0.5, width=0.99, y=0, height=0.93, just=c(0.5,0),
        xscale=c(0,sum(col.widths)), yscale=c(0, 100)))
      grid.rect(x=cumsum(col.widths)[unlist(lapply(set.sizes, function(s) c(rep(1,s),0)))[-(num.runs+num.traits)]==1]-0.5,
        y=rep(0, num.runs), height=at.least.one.rev.ordered.percent, width=0.5, gp=gpar(fill="gray70", col="gray50"),
        just=c(0.5,0), default.units="native")
      grid.segments(0, 0, 1, 0)
      grid.yaxis(at=c(0, 25, 50, 75, 100), gp=gpar(cex=0.6))
    popViewport()
    grid.text(c("percent","with at least","one reversal"), x=seq(-0.066, -0.040, length.out=3), y=0.5, rot=90, gp=gpar(cex=0.70))
  popViewport()
popViewport()
dev.off()


# END
