library(grid)
source("lib-closed.r")
source("~/code/colorfunc.r")

the.args <- commandArgs()
rimage <- sub("--rimage=", "", the.args[grep("--rimage=", the.args)])
pdf.file <- sub("--pdf=", "", the.args[grep("--pdf=", the.args)])

stopifnot(length(rimage)==1)
stopifnot(length(pdf.file)==1)
stopifnot(file.exists(rimage))
load(rimage)

tbl <- as.data.frame(t(sapply(run.res$fixations, unlist)))
sel.coef <- (tbl$fit[-1]/tbl$fit[-nrow(tbl)]-1) * sim.mod.start$N

use.cols <- grep("params", names(tbl))
param.names <- sub("params.", "", names(tbl)[use.cols])
x <- tbl[,"step"]
m <- run.res$mod.final

line.colors <- c("gray40", "purple1")
line.types <- c("solid", "11")
line.weights <- c(2, 3)
palette(c("white", "dodgerblue3", "gray50", "gold", "firebrick"))

plot.selection <- function(i) {
   y <- tbl[,i]
   y.rng <- expand.range2(y, amount=0.05)
   vp <- viewport(layout.pos.col=1, layout.pos.row=i, 
      xscale=c(0.5,length(y)+0.5), yscale=y.rng)
   pushViewport(vp)
   s <- sel.coef
   affect.this.parameter <- y[-1] != y[-length(y)] 
   sel.color <- rep(1, length(s))
   if (sum(affect.this.parameter) > 1) {
      sel.color[affect.this.parameter] <- apply(sapply(c(-2, 2, 5), 
         function(thresh) s[affect.this.parameter] > thresh), 1, sum)+2
   } else if (sum(affect.this.parameter) == 1) {
      sel.color[affect.this.parameter] <- sum(sapply(c(-2, 2, 5), 
         function(thresh) s[affect.this.parameter] > thresh))+2
   }
   sel.color <- c(sel.color, 1)

   # determine if we overshot
   overshot <- sapply(which(affect.this.parameter), function(j) {
      adj <- as.list(tbl[j,use.cols])
      names(adj) <- param.names
      q <- seq(adj[[i]], tbl[j+1,i], length.out=100)
      length(q) != which.max(sapply(q, function(x) {
         adj[[i]] <- x
         abs.fitness(do.call("alter.model", c(list(m),adj)))
      }))
   })
   # handle the case where there are no overshoots?
   if (length(overshot) == 0) overshot <- numeric(0)

   this.line.colors <- rep(line.colors[1], length(y))
   this.line.colors[-length(y)][affect.this.parameter][overshot] <- line.colors[2]
   this.line.types <- rep(line.types[1], length(y))
   this.line.types[-length(y)][affect.this.parameter][overshot] <- line.types[2]
   this.line.weights <- rep(line.weights[1], length(y))
   this.line.weights[-length(y)][affect.this.parameter][overshot] <- line.weights[2]
   
   a2col <- "gray30"
#   axis(3, col=a2col, col.ticks=a2col, col.axis=a2col)
#   mtext("time (mutations)", line=1.8, side=3, cex=0.75, col=a2col)
#   axis(1, at=seq(0,sim.mod.start$evolve$steps,length.out=6), labels=signif(seq(0, nrow(tbl), length.out=6), digits=1))
#   grid.segments(x2, tbl[,i], tbl$step, y.rng[2], gp=gpar(col="gray"))
   grid.rect(x=(1:length(y))-0.5, y.rng[2], width=1, height=y.rng[2]-y.rng[1], 
      gp=gpar(fill=c("gray85", "white"), col=c("gray85", "white")), 
      just=c("left", "top"), default.units="native")
   grid.rect(gp=gpar(col="gray50"))
#   rug(tbl$step, side=3, col=a2col)
   at <- pretty(y.rng, n=4)
   at <- at[y.rng[1] <= at & at <= y.rng[2]]
   grid.yaxis(at=at, gp=gpar(cex=0.7, col="gray30"))
   if (i == num.plots) {
      grid.xaxis(at=pretty(c(0, length(y)), n=4), gp=gpar(cex=0.8))
      grid.text("time (substitutions)", x=0.5, y=-0.4)
   }
   grid.segments(1:(length(y)-1), y[-length(y)], 2:length(y), y[-1], 
      default.units="native", gp=gpar(col=this.line.colors, lwd=this.line.weights, lty=this.line.types))
   grid.points(1:length(y), y, 
      size=unit(pt.sizes[(sel.color != 1 & sel.color != 3)+1], "char"),
      gp=gpar(fill=sel.color), pch=21)
   grid.text(param.expr[i], x=-0.2, y=0.5, gp=gpar(cex=1.2))
   popViewport()
}

pt.sizes <- c(0.5, 0.8)
num.plots <- length(sim.mod.start$evolve$params)
vp.plotset <- viewport(x=0.15, y=0.87, width=0.55, height=0.78, just=c("left", "top"),
   layout=grid.layout(nrow=num.plots, ncol=1))

# hard code in labels
stopifnot(sum(param.names != c("Y.a", "Y.b", "Z.a", "Z.b", "k"))==0)
param.expr <- c(expression(Y[alpha]), expression(Y[beta]), 
   expression(Z[alpha]), expression(Z[beta]), expression(k))

pdf(file=pdf.file, height=8, width=5.0)

# plot the main plots
pushViewport(vp.plotset)
trash <- lapply(1:num.plots, plot.selection)
popViewport()

# plot the substitution-mutation mapping at top
main.top <- 0.87
pushViewport(viewport(x=0.15, y=main.top, width=0.55, height=0.05, 
   xscale=c(0.5, nrow(tbl)+0.5), yscale=0:1, just=c("left", "bottom")))
x2 <- tbl$step/sim.mod.start$evolve$steps*(nrow(tbl)+1)
grid.segments(1:nrow(tbl), 0, x2, 1, default.units="native", gp=gpar(col="gray60"))
at <- pretty(c(0, sim.mod.start$evolve$steps), n=4)
grid.xaxis(at=at/sim.mod.start$evolve$steps*(nrow(tbl)+1), label=at, 
   main=FALSE, gp=gpar(cex=0.8))
grid.text("time (mutations)", x=0.5, y=2.20)
popViewport()

# plot the legend on the right
pushViewport(viewport(x=0.99, y=main.top, width=0.28, height=0.5, just=c("right", "top")))
ly <- seq(0.98, 0.3, length=6)
grid.segments(0.02, ly, 0.18, ly, 
   gp=gpar(col=line.colors[c(rep(1, 5), 2)], lwd=line.weights[c(rep(1, 5), 2)], lty=line.types[c(rep(1, 5), 2)]))
grid.points(x=rep(0.1, 5), y=ly[-6], pch=21, size=unit(pt.sizes[c(2,2,1,2,1)], "char"),
   gp=gpar(fill=c(5:1)))
grid.text(c(
      "beneficial", 
      "weakly beneficial",
      "effectively neutral",
      "deleterious",
      "mutation affects",
      "overshooting"), x=rep(0.2, 6), y=ly, just=c("left", "center"),
   gp=gpar(cex=0.87))
grid.text(c(
      "5 < Ns", 
      "2 < Ns < 5",
      "-2 < Ns < 2",
      "Ns < -2",
      "other parameter"), x=rep(0.2, 5), y=ly[-6]-0.05, just=c("left", "center"),
   gp=gpar(cex=0.87))
popViewport()

dev.off()


#par(mar=c(4,4,4,1), mgp=c(1.8,0.75,0), cex.lab=1.1, cex.axis=0.8)
#layout.m <- rep(0, plot.cols*plot.rows)
#layout.m[1:num.plots] <- 1:num.plots
#layout(matrix(layout.m, nrow=plot.rows, ncol=plot.cols, byrow=TRUE))
#trash <- lapply(1:length(param.names), function(i) {
#})
#plot(c(0,1), c(1,length(palette)), type="n", xlab="", ylab="", axes=FALSE)
#legend("topleft", pch=21, col=line.colors[c(rep(1,5),2)], pt.bg=palette(), pt.cex=c(0.8, 1.2, 0.8, 1.2, 1.2, NA), 
#   lwd=c(rep(NA, 5), 1),
#   legend=c(
#      "mutation affects other parameter",
#      "deleterious, Ns < -2",
#      "effectively neutral, -2 < Ns < 2",
#      "weakly beneficial, 2 < Ns < 5",
#      "beneficial, 5 < Ns", 
#      "overshooting"))
#   
#dev.off()

# END
