library(grid)
source("../src/r/multidimevol-lib.r")

# Require the argument --run=<runname> on the command line
the.args <- commandArgs()
run.name <- sub("--run=", "", the.args[grep("--run=", the.args)])
stopifnot(length(run.name) == 1)

# run.name for the submitted analyses:
# run.name <- "2011_10_20"

tol <- 1e-8

#param order in bitmask going from left to right:
#  k, Z.b, Z.a, Y.b, Y.a
ep.bitmask <- 31

false.density <- function(d) 24*d - 48*d^2 

# the fitness of a signal of duration, d
partial.fitness.sig <- function(d, mod) {
  y.level.max <- mod$Y$b / mod$Y$a * (1 - exp((-1) * mod$Y$a * d)) 
  cost.Y <- mod$Y$n * mod$Y$b / mod$Y$a * d
  if (y.level.max > mod$k) {
    delay <- (-1/mod$Y$a) * log(1 - mod$k * mod$Y$a / mod$Y$b)
    benefit <- mod$delta * mod$Z$b / mod$Z$a^2 * (mod$Z$a * (d - delay) - 1 + exp(-mod$Z$a*(d - delay)))
    cost.Z <- mod$Z$n * mod$Z$b / mod$Z$a * (d - delay)
  } else {
    benefit <- 0
    cost.Z <- 0
  }
  
  true.density <- 1/4 * mod$c2
  # compute the density of false signals
  if (d >= 1/2) {
    false.density <- 0
  } else {
    false.density <- mod$c1 * false.density(d)
  }

  true.density * (benefit - cost.Y - cost.Z) + false.density * ((-1) * (cost.Y + cost.Z))
}

m.start <- list(Y=list(a=2.05, b=0.95, n=0.4), Z=list(a=0.81, b=41.2, n=2), k=0.06, c1=0.05, c2=0.15, 
  delta=5, fit=list(offset=0, scale=1), timeresolution=0.01)
write.config.file(project.path(paste("analyses/configs/two_environments-start-", run.name, ".rconf", sep="")), m.start)

ev.em.low.noise <- evolve(100, m.start, tolerance=tol, max.step=0.10, evolvable.param.bitmask=ep.bitmask)
m.low.noise.optimum <- update.model.from.mx(m.start, ev.em.low.noise)
write.config.file(project.path(paste("analyses/configs/two_environments-low_noise_optimum-", run.name, ".rconf", sep="")), m.low.noise.optimum)

m.high.noise.optimum <- m.low.noise.optimum
m.high.noise.optimum$c1 <- 0.30
ev.em.high.noise <- evolve(100, m.high.noise.optimum, tolerance=tol, max.step=0.10, evolvable.param.bitmask=ep.bitmask)
m.high.noise.optimum <- update.model.from.mx(m.high.noise.optimum, ev.em.high.noise)
write.config.file(project.path(paste("analyses/configs/two_environments-high_noise_optimum-", run.name, ".rconf", sep="")), m.high.noise.optimum)

# only show the first 99.5% of the fitness increase
ev.em.high.noise <- ev.em.high.noise[
  ev.em.high.noise[,"fitness"] < (function(x) (x[2]-x[1])*0.995+x[1])(range(ev.em.high.noise[,"fitness"])),]

m.low.noise.opt.high.env <- m.low.noise.optimum
m.low.noise.opt.high.env$c1 <- m.high.noise.optimum$c1
m.high.noise.opt.low.env <- m.high.noise.optimum
m.high.noise.opt.low.env$c1 <- m.low.noise.optimum$c1

x <- seq(0, 4, 0.01)
low.opt.low.env <- sapply(x, partial.fitness.sig, m.low.noise.optimum)
low.opt.high.env <- sapply(x, partial.fitness.sig, m.low.noise.opt.high.env)
high.opt.high.env <- sapply(x, partial.fitness.sig, m.high.noise.optimum)
high.opt.low.env <- sapply(x, partial.fitness.sig, m.high.noise.opt.low.env)

fitnesses <- c(
  abs.fitness(m.low.noise.optimum), 
  abs.fitness(m.high.noise.opt.low.env),
  abs.fitness(m.low.noise.opt.high.env), 
  abs.fitness(m.high.noise.optimum)
)

# I also plot the delay, tau, over time during the evolution
ev.em.high.noise <- cbind(ev.em.high.noise, tau=with(as.data.frame(ev.em.high.noise), -1/Y.a * log(1-k*Y.a/Y.b)))


density.plot.max <- 1
pretty.range <- function(x) {
  xmin <- min(x)
  xmin.try <- signif(xmin * (1-1.4^(0:10)/100), digits=1)
  xmin <- xmin.try[which(xmin >= xmin.try)[1]]

  xmax <- max(x)
  xmax.try <- signif(xmax * (1+1.4^(0:10)/100), digits=1)
  xmax <- xmax.try[which(xmax <= xmax.try)[1]]
  return(c(xmin,xmax))
}

plot.height <- 3.75
plot.widths <- c(0.1, 0.3, 0.05, 0.3, 0.21, 0.1, 0.1, 0.05)
plot.heights <- c(0.07,0.3,0.13,0.14,0.04)
palette(c("#CD5B45", "gray30", rgb(0.7, 0.7, 0.7, 0.5), rgb(205, 91, 69, 128, maxColorValue=255), 
  "gray15", rgb(0.4, 0.4, 0.4, 0.5)))
pdf(file=project.path(paste("analyses/plots/two_environments-", run.name, ".pdf", sep="")), 
  height=plot.height, width=plot.height*sum(plot.widths)/sum(plot.heights))
pushViewport(viewport(layout=grid.layout(nrow=length(plot.heights), ncol=length(plot.widths), 
  widths=plot.widths, heights=plot.heights)))
y.rng <- expand.range2(range(c(low.opt.low.env, low.opt.high.env, high.opt.low.env, high.opt.high.env)))
# low noise environment
pushViewport(viewport(layout.pos.row=2, layout.pos.col=2,
  xscale=range(x), yscale=c(0,density.plot.max)))
grid.polygon(x[x <= 1/2], m.low.noise.optimum$c1 * false.density(x[x <= 1/2]), 
  default.units="native", gp=gpar(fill=3, col=2))
grid.polygon(c(x[1],x,x[length(x)]), m.low.noise.optimum$c2 * c(0, rep(1/4, length(x)), 0), 
  default.units="native", gp=gpar(fill=6, col=5))
grid.text("signal duration", x=0.5, y=-0.20, gp=gpar(cex=0.85))
grid.text("partial fitness", x=-0.25, y=0.5, rot=90, gp=gpar(cex=0.85))
popViewport()
pushViewport(viewport(layout.pos.row=2, layout.pos.col=2))
legend.y <- seq(0.45, 0.15, length.out=4)
legend.x <- 0.45
grid.text(c("low-noise adapted", "high-noise adapted", "false signal distribution", "true signal distribution"),
  x=legend.x, y=legend.y, gp=gpar(cex=0.65), just=c(0,0.5))
grid.rect(x=legend.x-0.06, y=legend.y[3:4], height=0.06, width=0.06, gp=gpar(fill=c(3,6), col=c(2,5)))
grid.segments(legend.x-0.09, legend.y[1:2], legend.x-0.03, legend.y[1:2], gp=gpar(col=c(1,2), lwd=2))
popViewport()
# low noise partial fitness
pushViewport(viewport(layout.pos.row=2, layout.pos.col=2, xscale=range(x), yscale=y.rng))
grid.lines(x, low.opt.low.env, default.units="native", gp=gpar(col=1, lwd=2))
grid.lines(x, high.opt.low.env, default.units="native", gp=gpar(col=2, lwd=2))
grid.xaxis(gp=gpar(cex=0.7, lineheight=0.8))
grid.yaxis(gp=gpar(cex=0.7, lineheight=0.8))
popViewport()
# high noise environment
pushViewport(viewport(layout.pos.row=2, layout.pos.col=4,
  xscale=range(x), yscale=c(0,density.plot.max)))
grid.polygon(x[x <= 1/2], m.high.noise.optimum$c1 * false.density(x[x <= 1/2]), 
  default.units="native", gp=gpar(fill=3, col=2))
grid.polygon(c(x[1],x,x[length(x)]), m.high.noise.optimum$c2 * c(0, rep(1/4, length(x)), 0), 
  default.units="native", gp=gpar(fill=6, col=5))
grid.yaxis(main=FALSE, gp=gpar(cex=0.7, lineheight=0.8))
grid.text("signal duration", x=0.5, y=-0.20, gp=gpar(cex=0.85))
grid.text("signal density", x=1.25, y=0.5, rot=-90, gp=gpar(cex=0.85))
popViewport()
# high noise partial fitness
pushViewport(viewport(layout.pos.row=2, layout.pos.col=4, xscale=range(x), yscale=y.rng))
grid.lines(x, low.opt.high.env, default.units="native", gp=gpar(col=1, lwd=2))
grid.lines(x, high.opt.high.env, default.units="native", gp=gpar(col=2, lwd=2))
grid.xaxis(gp=gpar(cex=0.7, lineheight=0.8))
popViewport()
# fitness plot
pushViewport(viewport(layout.pos.row=2, layout.pos.col=6,
  xscale=c(0.5,2.5), yscale=expand.range2(fitnesses, amount=0.2)))
grid.text("total fitness", x=-0.6, y=0.5, rot=90, gp=gpar(cex=0.85))
grid.yaxis(gp=gpar(cex=0.7, lineheight=0.8))
grid.points(1:2, fitnesses[1:2], default.units="native", pch=20, gp=gpar(col=1:2))
grid.segments(unit(1, "native"), 0.94, unit(1, "native"), 1.03, gp=gpar(col=1))
grid.segments(unit(2, "native"), 0.94, unit(2, "native"), 1.12, gp=gpar(col=2))
grid.segments(unit(1, "native"), 1.03, 1, 1.03, gp=gpar(col=1))
grid.segments(unit(2, "native"), 1.12, 1, 1.12, gp=gpar(col=2))
grid.rect(0.7, 1.07, width=1.5, height=0.04, gp=gpar(col=rgb(1,1,1,0.5), fill=rgb(1,1,1,0.5)))
grid.text("low-noise adapted", 0.7, 1.08, gp=gpar(col=1, cex=0.6))
grid.rect()
popViewport()
pushViewport(viewport(layout.pos.row=3, layout.pos.col=6))
grid.text(c("low noise", "envir."), x=0.5, y=c(0.85, 0.65), gp=gpar(cex=0.63))
popViewport()
pushViewport(viewport(layout.pos.row=3, layout.pos.col=7))
grid.text(c("high noise", "envir."), x=0.5, y=c(0.85, 0.65), gp=gpar(cex=0.63))
popViewport()
pushViewport(viewport(layout.pos.row=2, layout.pos.col=7,
  xscale=c(0.5,2.5), yscale=expand.range2(fitnesses, amount=0.2)))
grid.points(1:2, fitnesses[3:4], default.units="native", pch=20, gp=gpar(col=1:2))
grid.text("high-noise adapted", 0.4, 1.17, gp=gpar(col=2, cex=0.6))
grid.segments(unit(1, "native"), 0.94, unit(1, "native"), 1.03, gp=gpar(col=1))
grid.segments(unit(2, "native"), 0.94, unit(2, "native"), 1.12, gp=gpar(col=2))
grid.segments(0, 1.03, unit(1, "native"), 1.03, gp=gpar(col=1))
grid.segments(0, 1.12, unit(2, "native"), 1.12, gp=gpar(col=2))
grid.rect()
popViewport()
# panel labels
pushViewport(viewport(layout.pos.row=1, layout.pos.col=1))
grid.text("A", x=0.3)
popViewport()
pushViewport(viewport(layout.pos.row=1, layout.pos.col=5))
grid.text("B", x=0.6)
popViewport()
pushViewport(viewport(layout.pos.row=3, layout.pos.col=1))
grid.text("C", x=0.3, y=0.30)
popViewport()
pushViewport(viewport(layout.pos.row=3, layout.pos.col=6))
grid.text("D", x=0.5, y=0.30)
popViewport()
# lower row of param-over-time plots
pushViewport(viewport(layout.pos.row=4, layout.pos.col=c(1,length(plot.widths))))
param.names <- c("Y.a", "Y.b", "Z.a", "Z.b", "k")
pushViewport(viewport(x=0.99, width=0.94, 
  layout=grid.layout(nrow=1,ncol=length(param.names)+3, widths=c(rep(1,6), 0.32, 1)), just=c(1,0.5)))
time.subset <- floor(seq(1, nrow(ev.em.high.noise), length.out=100))
stopifnot(sum(colnames(ev.em.high.noise)[3:7] != param.names) == 0)
param.expr <- c(expression(alpha[y]), expression(beta[y]),
   expression(alpha[z]), expression(beta[z]), "k")
# standarize the scale of the first three parameters
param.ranges <- apply(ev.em.high.noise, 2, function(x) pretty.range(expand.range2(range(x))))
#use.common <- c("Y.a", "Y.b", "Z.a")
#common.range <- range(param.ranges[,use.common])
#param.ranges[,use.common] <- common.range
# manual param ranges:
param.ranges[,"Y.a"] <- c(0.6,2)
param.ranges[,"Y.b"] <- c(0.6,2)
param.ranges[,"Z.a"] <- c(0.6,2)
param.ranges[,"Z.b"] <- c(40,50)
param.ranges[,"k"] <- c(0.03,0.4)
label.offset <- c(-0.05, -0.050, -0.035, -0.035, 0.02)
for (i in 1:length(param.names)) {
  pushViewport(viewport(layout.pos.row=1, layout.pos.col=i, 
    xscale=range(time.subset), yscale=param.ranges[,param.names[i]]))
  xl <- seq(0,1,len=5)
  grid.segments(xl, 0, xl, 1, gp=gpar(col="gray80"))
  grid.segments(0, xl, 1, xl, gp=gpar(col="gray80"))
  grid.lines(time.subset, ev.em.high.noise[time.subset,param.names[i]], 
    default.units="native", gp=gpar(col=2, lwd=2))
#  grid.segments(0, 0, 0:1, c(1.1,0))
  grid.segments(0, 0, 0, 1.1)
  grid.text(param.ranges[,param.names[i]], x=rep(-0.05,2), y=c(0,1), 
    just=c(1,0.5), gp=gpar(cex=0.6, col="gray30"))
  grid.segments(0, 0, -0.03, 0)
  grid.segments(0, 1, -0.03, 1)
  grid.text(param.expr[i], x=0.05, y=1.0+label.offset[i], just=c(0,0))
  popViewport()
}
# fitness plot
pushViewport(viewport(layout.pos.row=1, layout.pos.col=length(param.names)+1, 
  xscale=range(time.subset), yscale=param.ranges[,"fitness"]))
xl <- seq(0,1,len=5)
grid.segments(xl, 0, xl, 1, gp=gpar(col="gray80"))
grid.segments(0, xl, 1, xl, gp=gpar(col="gray80"))
grid.lines(time.subset, ev.em.high.noise[time.subset,"fitness"], 
  default.units="native", gp=gpar(col=1, lwd=2))
grid.segments(0, 0, 0, 1.1)
grid.text("fitness", x=0.03, y=1.15, just=c(1,1), rot=90, gp=gpar(cex=0.75, col=1))
grid.text(param.ranges[,"fitness"], x=rep(-0.05,2), y=c(0,1), 
  just=c(1,0.5), gp=gpar(col="gray30", cex=0.6))
grid.segments(0, 0, -0.03, 0)
grid.segments(0, 1, -0.03, 1)
popViewport()
# tao plot
pushViewport(viewport(layout.pos.row=1, layout.pos.col=length(param.names)+3, 
  xscale=range(time.subset), yscale=expand.range2(range(ev.em.high.noise[,"tau"]))))
xl <- seq(0,1,len=5)
grid.segments(xl, 0, xl, 1, gp=gpar(col="gray80"))
grid.segments(0, xl, 1, xl, gp=gpar(col="gray80"))
grid.lines(time.subset, ev.em.high.noise[time.subset,"tau"], 
  default.units="native", gp=gpar(col=2, lwd=2))
grid.segments(0, 0, 0, 1.1)
grid.text(expression(tau), x=0.05, y=1.07, just=c(0,0.5))
grid.text(param.ranges[,"tau"], x=rep(-0.05,2), y=c(0,1), 
  just=c(1,0.5), gp=gpar(cex=0.6, col="gray30"))
grid.segments(0, 0, -0.03, 0)
grid.segments(0, 1, -0.03, 1)
popViewport()
popViewport()
popViewport()
pushViewport(viewport(layout.pos.row=5, layout.pos.col=c(1,length(plot.widths))))
grid.segments(0.3, 0.3, 0.7, 0.3, arrow=arrow(type="closed", length=unit(0.25, "npc")), gp=gpar(fill="black", lwd=1.3))
grid.rect(0.5, 0.35, width=0.06, height=0.7, gp=gpar(col="white", fill="white"))
grid.text("time", y=0.4, gp=gpar(cex=0.85))
popViewport()
# cleanup
popViewport()
dev.off()


# END


