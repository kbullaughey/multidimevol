# sample from a hypercube
sample.choices <- function(n) {
  data.frame(
    Y.a=runif(n, min=0.5, max=3),
    Y.b=runif(n, min=0.25, max=1.5),
    Z.a=runif(n, min=0.25, max=1.5),
    Z.b=runif(n, min=30, max=60),
    k=runif(n, min=0.05, max=0.5))
}
choices.mx <- sample.choices(10000)

m.start <- load.config(file="../configs/high_noise_optimum.rconf", defaults=default.model)

