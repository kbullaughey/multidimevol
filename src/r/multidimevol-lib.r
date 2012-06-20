# The project must be contained in a directory labeled 'multidimevol' otherwise
# scripts won't be able to find the project root
project.root <- gsub("/multidimevol/.*$", "/multidimevol", getwd())
project.path <- function(relpath) paste(project.root, relpath, sep="/")

dyn.load(project.path("src/c/reg_xyz.so"))

zip <- function(x, y, func) sapply(1:length(x), function(i) func(x[[i]], y[[i]]))
reduce.vector <- function(x, len=NULL) {
   if (is.null(len)) {
      return(x)
   }
   x[round(seq(1, length(x), length.out=len))]
}

fitness <- function(m) {
   params <- as.double(c(
      m$Y$a, m$Y$b, m$Y$n,
      m$Z$a, m$Z$b, m$Z$n,
      m$k, m$c1, m$c2, m$delta))
   fitness <- 0
   error <- 0
   res <- .C("c_fitness", m_in=params, m_len=as.integer(length(params)), 
      error=as.integer(error), fit=as.double(fitness))
   if (res$error!=0) {
      cat("ERROR:", res$error, "\n")
      return(NULL)
   }
   return(res$fit)
}

abs.fitness <- function(m) {
   w <- fitness(m)
   pmax(0, w*m$fit$scale + m$fit$offset)
}


# points is a matrix with 5 columns corresponding to the parameters and rows 
#  corresponding to sets of points
vector.field <- function(points, m) {
   stopifnot(ncol(points) == 5 && nrow(points) >= 1)
   param.mx <- as.double(t(points))
   params <- as.double(c(
      m$Y$a, m$Y$b, m$Y$n,
      m$Z$a, m$Z$b, m$Z$n,
      m$k, m$c1, m$c2, m$delta))
   vector.field <- double(length(param.mx))
   error <- 0

   res <- .C("c_gradient_field", m_in=params, m_len=as.integer(length(params)), 
      param_matrix=param.mx, pm_len=as.integer(length(param.mx)),
      vector_field=vector.field, vf_len=as.integer(length(vector.field)),
      error=as.integer(error))
   if (res$error!=0) {
      cat("ERROR:", res$error, "\n")
      return(NULL)
   }
   mx <- matrix(res$vector_field, ncol=5, nrow=res$vf_len/5, byrow=TRUE)
   colnames(mx) <- c("Y.a", "Y.b", "Z.a", "Z.b", "k")
   return(mx)
}

evolve <- function(dur, m, tolerance=0.1, max.step=-1, evolvable.param.bitmask=NULL) {
   evolvable.params <- 5
   if (is.null(evolvable.param.bitmask)) {
      evolvable.param.bitmask <- 2^evolvable.params - 1
   }
   stopifnot(evolvable.param.bitmask < 2^evolvable.params)
   param.evol <- double((2+evolvable.params)*ceiling(dur/m$timeresolution+2))
   params <- as.double(c(
      m$Y$a, m$Y$b, m$Y$n,
      m$Z$a, m$Z$b, m$Z$n,
      m$k, m$c1, m$c2, m$delta))

   error <- 0
   out.rows <- 0
   res <- .C("c_evolution_matrix", m_in=params, m_len=as.integer(length(params)), 
      pm=param.evol, pm_len=as.integer(length(param.evol)), pm_rows_out=as.integer(out.rows), 
      max_duration=as.double(dur), time_resolution=as.double(m$timeresolution), tolerance=as.double(tolerance), 
      max_step=as.double(max.step), error=as.integer(error), evolvable_param_bitmask=as.integer(evolvable.param.bitmask))
   mx <- matrix(res$pm[1:((2+evolvable.params)*res$pm_rows_out)], 
      nrow=res$pm_rows_out, ncol=2+evolvable.params, byrow=TRUE)
   if (res$error!=0) {
      cat("ERROR:", res$error, "\n")
      return(NULL)
   }
   colnames(mx) <- c("time", "fitness", "Y.a", "Y.b", "Z.a", "Z.b", "k")
   return(mx)
}

model.contains <- function(m, addresses) 
   sapply(model.lookup(m, addresses), function(x) !is.null(x))

model.address <- function(addresses, single=FALSE) {
   a <- strsplit(addresses, split="\\.")
   if (single) {
      return(a[[1]])
   } 
   return(a)
}
model.update <- function(m, address, value) {
   cmd <- paste(paste(c("m", address), collapse="$"), "<-", deparse(value))
   eval(parse(text=cmd))
   return(m)
}
model.lookup <- function(m, addresses, single=FALSE) {
   if (single) {
      stopifnot(length(addresses)==1)
      a <- addresses
      s <- model.address(a, single=TRUE)
      res <- eval(parse(text=paste("m", paste("[[\"", s, "\"]]", sep="", collapse=""), sep="")))
   } else {
      res <- lapply(addresses, function(a) {
         s <- model.address(a, single=TRUE)
         eval(parse(text=paste("m", paste("[[\"", s, "\"]]", sep="", collapse=""), sep="")))
      })
   }
   names(res) <- addresses
   return(res)
}

alter.model.list <- function(m, args) {
   addresses <- model.address(names(args))
   for (i in 1:length(args)) {
      m <- model.update(m, addresses[[i]], args[[i]])
   }
   return(m)
}

alter.model <- function(m, ...) {
   args <- list(...)
   return(alter.model.list(m, args))
}

update.model.from.mx <- function(m, mx, row=nrow(mx)) {
   to.update <- as.list(mx[row,])
   alter.model.list(m, to.update[model.contains(m, names(to.update))])
}

expanded.range <- function(y) range(c(y, y[1]*c(0.9,1.1), y[length(y)][c(0.9,1.1)]))
expand.range2 <- function(y, amount=0.1) range(c(y, min(y)*(1 + amount*c(-1,1)), max(y)*(1 + amount*c(-1,1))))

collate.par.vs.time <- function(runs, max.time=Inf) {
   lapply(runs, function(ev) {
      ev <- ev[ev[,"time"] < max.time,]
      ev[floor(seq(1,nrow(ev),length.out=200)),]
   })
}

plot.par.vs.time <- function(evc, prefix=NULL, max.time=Inf, ranges=NULL, gray.lines=TRUE, 
      interior.points=NULL) {
   par(mar=c(4,4,1,1), mgp=c(2,0.8,0), cex.lab=2.3, cex.axis=1.5)
   palette(c("firebrick", rep("purple3", ncol(evc[[1]])-2)))
   trash <- lapply(1:length(evc), function(r) {
      if (!is.null(prefix)) {
         png(file=paste(prefix, "-", r, ".png", sep=""), height=700, width=950)
      }
      layout(matrix(1:6, nrow=2))
      ev <- evc[[r]]
      lapply(2:ncol(ev), function(i) {
         lab <- colnames(ev)[i]
         y <- ev[,i]
         if (!is.null(ranges) && !is.null(ranges[[lab]])) {
            y.rng <- ranges[[lab]]
         } else {
            y.rng <- expanded.range(y)
         }
         plot(range(ev[,1]), y.rng, type="n", xlab="time", ylab=lab)
         if (gray.lines) {
            trash2 <- lapply(1:length(evc), function(j) {
               lines(reduce.vector(evc[[j]][,1], len=interior.points), 
                  reduce.vector(evc[[j]][,i], len=interior.points), col="gray", lwd=2)
            })
         }
         lines(reduce.vector(ev[,1], len=interior.points), 
            reduce.vector(y, len=interior.points), col=i-1, lwd=5)
      })
      if (!is.null(prefix)) dev.off()
   })
}

default.model <- function(...) {
   alter <- list(...)
   args <- list(list(),
      Y.a=2.5, Y.b=1.0, Y.n=0.2,
      Z.a=2.5, Z.b=1.0, Z.n=0.2,
      k=0.2, c1=0.3, c2=0.05, delta=7,
      fit.offset=2, fit.scale=1, timeresolution=0.01)
   m <- do.call(alter.model, args)
   if (length(alter) > 0) {
      m <- do.call(alter.model, c(list(m),alter))
   }
   return(m)
}

default.simulation <- function(...) {
   alter <- list(...)
   args <- list(list(),
      N=1000, evolve.steps=10000, evolve.params=c("Y.b", "Z.b", "k"), 
      mut.sd=1/50, mut.components=1, evolve.boundaries=NULL, mut.num=1,
      mut.func=function(x,delt) x+delt)
   s <- do.call(alter.model, args)
   if (length(alter) > 0) {
      s <- do.call(alter.model, c(list(s),alter))
   }
   return(s)
}

# delist always returns a named list
delist <- function(x) {
   if (is.list(x)) {
      hn <- names(x)[1]
      if (length(x) > 1) {
         h <- x[[1]]
         t <- x[-1]
         dt <- delist(t)
         if (is.list(h)) {
            dh <- delist(h)
            names(dh) <- paste(hn, names(dh), sep=".")
            res <- c(dh, dt)
         } else {
            res <- c(list(h), dt)
            names(res)[1] <- hn
         }
      } else {
         x <- x[[1]]
         if (is.list(x)) {
            res <- delist(x)
            names(res) <- paste(hn, names(res), sep=".")
         } else {
            res <- list(x)
            names(res) <- hn
         }
      }
   } else {
      res <- list(x)
      names(res) <- hn
   }
   return(res)
}

write.config.file <- function(file, m) {
   x <- delist(m)
   cat(file=file, "# Created: ", date(), "\n", sep="")
   cat(file=file, paste(names(x), x, sep=": "), sep="\n", append=TRUE)
}

# the argument, defaults, is a function to build the default configuration
#  this function is used both for simulation configurations and model configurations
load.config <- function(file, defaults=default.model) {
   lines <- scan(file, what="", sep="\n", comment.char="#")
   conf.mx <- matrix(unlist(strsplit(lines, split=": ")), ncol=2, byrow=TRUE)
   args <- lapply(conf.mx[,2], function(x) eval(parse(text=x)))
   # make sure no parameters are invalid
   b <- grep("^evolve.boundaries", conf.mx[,1])
   stopifnot(sum(!(conf.mx[-b,1] %in% names(delist(defaults())))) == 0)
   names(args) <- conf.mx[,1]
   m <- do.call(defaults, args)
   return(m)
}

# fixation probability
prob.fix <- function(mut, wt, N) {
   stopifnot(length(mut) == length(wt) || length(wt) == 1)
   stopifnot(sum(wt==0)==0)
   s <- mut/wt - 1
   pfix <- (1 - exp(-s))/(1 - exp(-N*s))
   pfix[wt == mut] <- 1/N
   pfix[mut == 0] <- 0
   return(pfix)
}

summarize.simulation <- function(to.mutate, p, fit, i) {
  cat(paste(zip(p, to.mutate, function(x, y) paste(signif(x, digits=5), c(" ", "*")[y+1], sep="")), sep="\t"), i, fit, sep="\t")
  cat("\n")
}

popsim <- function(mod.sim, mod.reg) {
   r <- length(mod.sim$evolve$params)
   stopifnot(mod.sim$mut$num <= r && mod.sim$mut$num > 0)
   # check the boundary specification if there is one
   boundaries <- model.lookup(mod.sim, paste("evolve.boundaries.", mod.sim$evolve$params, sep=""))
   names(boundaries) <- sub("evolve.boundaries.", "", names(boundaries))
   no.boundary <- unique(c(grep("empty$", names(boundaries))))
   no.boundary <- sapply(1:length(boundaries), function(x) x %in% no.boundary)
   sel <- sapply(boundaries, is.null)
   boundaries[sel & !no.boundary] <- rep(list(c(0,Inf)), sum(sel & !no.boundary))
   boundaries[sel & no.boundary] <- rep(list(c(-Inf,Inf)), sum(sel & no.boundary))
   rm(sel, no.boundary)
   stopifnot(length(mod.sim$mut$sd)==length(mod.sim$mut$components))
   stopifnot(sum(model.contains(mod.reg, mod.sim$evolve$params)) == r)
   wt <- mod.reg
   wt.fit <- abs.fitness(wt)
   fixations <- list()
   pwt <- unlist(model.lookup(wt, mod.sim$evolve$params))
   fixations[[1]] <- c(params=pwt, step=0, fit=wt.fit)
   cat(paste(pwt, sep="\t"), 0, wt.fit, sep="\t")
   cat("\n")
   for (i in 1:mod.sim$evolve$steps) {
      x <- as.numeric(pwt)
      to.mutate <- rep(0, r)
      to.mutate[sample(1:r, mod.sim$mut$num, replace=FALSE)] <- 1
      mut.sd <- mod.sim$mut$sd[sample(1:length(mod.sim$mut$components), r, prob=mod.sim$mut$components, replace=TRUE)]
      delta <- rnorm(r, mean=0, sd=mut.sd) * to.mutate
      # assume reflective boundaries, here I make the assumption that that 
      x2 <- mod.sim$mut$func(x,delta)
      flip <- zip(x2, boundaries, function(z, b) z < b[1] || b[2] < z)
      delta[flip] <- delta[flip] * (-1)
      x3 <- mod.sim$mut$func(x,delta)
      flip <- zip(x3, boundaries, function(z, b) z < b[1] || b[2] < z)
      if (sum(flip) > 0) { 
         cat("still out of range\n", file=stderr())
         next
      }

      mut <- wt
      x <- x3
      for (k in 1:r) {
         mut <- model.update(mut, model.address(mod.sim$evolve$params[k], single=TRUE), x[k])
      }
      mut.fit <- abs.fitness(mut)

      p <- prob.fix(mut.fit, wt.fit, mod.sim$N)
      if (runif(1) < p) {
         # accept the mutation as fixing
         wt <- mut
         wt.fit <- mut.fit
         pwt <- unlist(model.lookup(wt, mod.sim$evolve$params))
         fixations[[length(fixations)+1]] <- c(params=pwt, step=i, fit=wt.fit)
         summarize.simulation(to.mutate, pwt, wt.fit, i)
      }
   }
   return(list(fixations=fixations, mod.final=wt))
}

bf <- function(bool) bool+0
evolvable.traits <- function(m) c(m$Y$a, m$Y$b, m$Z$a, m$Z$b, m$k)

fitness.confirm <- function(r) {
  # Currently I haven't yet abstracted the ranges for the signal distributions,
  #   so I manually add these into the model. Later these could be incorporated
  # I also change the locations of c1 and c2 to match the manuscript notation.
  r$c <- list(F=r$c1, T=r$c2)
  r$t <- list(F=1/2, T=4)
  
  y.cost <- r$Y$b * r$Y$n * (r$c$F * r$t$F + r$c$T * r$t$T) / (2*r$Y$a)

  if (r$k * r$Y$a / r$Y$b > 1) {
    z.cost.false <- 0
    z.cost.true <- 0
    z.benefit <- 0
  } else {
    tau <- (-1)/r$Y$a * log(1 - r$k * r$Y$a / r$Y$b)
    z.cost.false <- r$c$F * bf(tau < r$t$F) * 
      (r$Z$b * r$Z$n) / (2 * r$Z$a * r$t$F^3) * (r$t$F + tau) * (r$t$F - tau)^3
    z.cost.true <- r$c$T * bf(tau < r$t$T) * 
      (r$Z$b * r$Z$n) / (2 * r$Z$a * r$t$T) * (r$t$T - tau)^2
    z.benefit <- r$c$T * bf(tau < r$t$T) * 
      (r$Z$b * r$delta) / (r$Z$a^3 * r$t$T) * 
      (1 - r$Z$a * r$t$T - exp(-r$Z$a * (r$t$T - tau)) + r$Z$a * tau + (1/2)*r$Z$a^2*(r$t$T - tau)^2)
  }

  penalty.low <- exp((1/25) * sum(evolvable.traits(r)))
  penalty.high <- exp((-1)*sum(log(evolvable.traits(r))))
  
  return((-1)*y.cost - z.cost.false - z.cost.true + z.benefit - penalty.low - penalty.high)
}

# END



