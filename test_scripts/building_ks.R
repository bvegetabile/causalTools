source('~/git/causalTools/causalTools.R')

n_obs <- 10000
X <- rnorm(n_obs)
ps <- pnorm(0.5*X)
TA <- rbinom(n_obs, 1, ps)
XT <- X[TA==1]
XC <- X[TA==0]

plot(XT, f1(XT))

my_cdf <- function(dat, pts){
  n_pts <- length(pts)
  results <- matrix(NA, nrow = length(pts), ncol=2)
  for(d in 1:n_pts){
    results[d, ] <- c(pts[d], mean(dat < pts[d]))
  }
  return(results)
}

R_cdf <- function(dat, pts){
  fn <- ecdf(dat)
  return(fn(pts))
}

test_pts <- rnorm(10000)
out1 <- my_cdf(X, test_pts)
out2 <- R_cdf(X, test_pts)
out3 <- c_ecdf(X, test_pts)
out4 <- c_ecdf2(X, test_pts)
all.equal(out1[,2], out2)
all.equal(out3[,2], out2)
all.equal(out4[,2], out2)

library(rbenchmark)
benchmark(my_cdf(X, test_pts), 
          R_cdf(X, test_pts),
          c_ecdf(X, test_pts),
          c_ecdf2(X, test_pts),
          bal_stats(X, TA),
          replications=10,
          order='relative')

ps

t_wts <- (TA / ps) / sum(TA / ps)
c_wts <- ((1-TA) / (1-ps)) / sum(((1-TA) / (1-ps)))
plot(wtd_ecdf(X, wts = t_wts))
lines(wtd_ecdf(XT, (1/ps)[TA==1]), col='blue')
lines(wtd_ecdf(X, wts = c_wts), col='red')
lines(wtd_ecdf(X, rep(1/n_obs, n_obs)))

t_fn <- wtd_ecdf(X, wts = t_wts)
c_fn <- wtd_ecdf(X, wts = c_wts)
t_fn(-10:10) - c_fn(-10:10)

benchmark(mean(ks_avg_test(X, TA, ps, 100)))

source('~/git/causalTools/causalTools.R')

n_obs <- 10000
X <- rbinom(n_obs,1,0.3)
ps <- pnorm(0.5*X)
TA <- rbinom(n_obs, 1, ps)
XT <- X[TA==1]
plot(X, ps)
mean(ks_avg_test(X, TA, rep(1/n_obs, n_obs), 23509))
t_wts <- (TA / ps) / sum(TA / ps)
c_wts <- ((1-TA) / (1-ps)) / sum(((1-TA) / (1-ps)))
plot(wtd_ecdf(X, wts = t_wts))
lines(wtd_ecdf(XT, (1/ps)[TA==1]), col='blue')

plot(wtd_ecdf(X, wts = rep(1/n_obs, n_obs)))
lines(wtd_ecdf(XT, rep(1/n_obs, n_obs)[TA==1]), col='blue')

