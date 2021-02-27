############################################################################################################################################################################################
# Main file for the simulation experiments for the P\'olya-Gamma sampler of a generalized logistic regression
# Dalla Valle, L., Leisen, F., Rossini, L. and Zhu, W. (2020) - "A Polya-Gamma sampler for a generalized logistic regression"
############################################################################################################################################################################################

rm(list = ls())

# Load the necessary library

library(MASS)
library(truncnorm)
library(Rcpp)
library(RcppArmadillo)

# Set the directory where save the results of the Gibbs sampler
# setwd("")
sourceCpp("Fabrizio.cpp")

##############################################################################################
##                            SLICE SAMPLER FUNCTION
##############################################################################################
sliceSample = function (n, logf, x.interval = c(0, 3), root.accuracy = 0.01) {
  # n is the number of points wanted
  # f is the target distribution
  # x.interval is the A,B range of x values possible.
  
  pts = vector("numeric", n)  # This vector will hold our points.
  x = runif(1, x.interval[1], x.interval[2]) # Take a random starting x value.
  for (i in 1:n) {
    pts[i] = x
    y = logf(x) - rexp(1,1)    # Take a random y value
    # Imagine a horizontal line across the distribution. Find intersections across that line.
    fshift = function (x) { logf(x) - y }
    roots = c()
    for (j in seq(x.interval[1] + root.accuracy, x.interval[2] - root.accuracy, by = root.accuracy)) {
      if(!is.na(fshift(j)) & !is.na(fshift(j + root.accuracy))){
        if ((fshift(j) < 0) != (fshift(j + root.accuracy) < 0)) {
          # Signs don't match, so we have a root.
          root = uniroot(fshift, c(j, j + root.accuracy))$root
          roots = c(roots, root)
        }
      }
    }
    # Include the endpoints of the interval.
    roots = c(x.interval[1], roots, x.interval[2])
    
    if(length(roots)==2){# No roots found!!!
      warning("no roots found")
      x = runif(1, x.interval[1], x.interval[2]) # Take a random starting x value.
    }else{
      # Divide intersections into line segments.
      segments = matrix(ncol = 2)
      for (j in 1:(length(roots) - 1)) {
        midpoint = (roots[j + 1] + roots[j]) / 2.0
        if (logf(midpoint) > y) {
          # Since this line segment is under the curve, add it to segments.
          segments = rbind(segments, c(roots[j], roots[j + 1]))
        }
      }
      # Uniformly sample next x from segments
      # Assign each segment a probability, then unif based on those probabilities.
      # This is a bit of a hack because segments first row is NA, NA
      # Yet, subsetting it runs the risk of reducing the matrix to a vector in special case.
      total = sum(sapply(2:nrow(segments), function (i) {
        segments[i, 2] - segments[i, 1]
      }))
      probs = sapply(2:nrow(segments), function (i) {
        (segments[i, 2] - segments[i, 1]) / total
      })
      # Assign probabilities to each line segment based on how long it is
      # Select a line segment by index (named seg)
      p = runif(1, 0, 1)
      selectSegment = function (x, i) {
        if (p < x) return(i)
        else return(selectSegment(x + probs[i + 1], i + 1))
      }
      seg = selectSegment(probs[1], 1)
      
      # Uniformly sample new x value
      x = runif(1, segments[seg + 1, 1], segments[seg + 1, 2])
    }
  }
  return(pts)
}

##############################################################################################
##                          LOG POSTERIOR OF P
# Full joint log posterior of p for the generalized logistic regression model:
# prior for p ~ Gamma(hyper.a, hyper.b)
##############################################################################################
log_condi_p<-function(p, z, xbeta, hyper.a = 1, hyper.b = 1){
  n      <- length(z)
  xbeta <- ifelse(xbeta>10, 10, xbeta)
  like   <- -n*log(beta(p, p)) + sum(p*(z-xbeta) - 2*p*log(1+exp(z-xbeta)))
  prior  <- dgamma(p, hyper.a, hyper.b, log = TRUE)
  like+prior
}

##############################################################################################
##                          LOG POSTERIOR OF P
# Full joint log posterior of p for the generalized logistic regression model given y
# prior for p ~ Gamma(hyper.a, hyper.b)
##############################################################################################
log_condi_p_y <- function(p, y, xbeta, hyper.a = 1, hyper.b = 1){
  n <- length(y)
  xbeta <- ifelse(xbeta>10, 10, xbeta)
  prob.p <- pbeta(expit(xbeta), p, p)
  like <- sum(dbinom(y, 1, prob.p, log=TRUE))
  prior  <- dgamma(p, hyper.a, hyper.b, log = TRUE)
  like + prior
}

##############################################################################################
##                          GIBBS SAMPLING 
##############################################################################################
Bayes.logistic<-function(y,X,p,beta0,
                         prior.mn=0, prior.sd=10, hyper.a, hyper.b,
                         n.samples,burn,n.slice=10){
  # MCMC code for the model:
  can.sd=0.4
  nc <- ncol(X)
  n <- length(y)
  #Initial values:
  xbeta <- X %*% beta0
  w <- rep(0,n)
  for (j in 1:n) {
    w[j] <- rpg(2*p, 0) 
  }
  low = rep(-Inf, n)
  low[y==1] = 0
  up = rep(Inf, n)
  up[y==0] = 0 # upper and lower boundry for truncated normal
  
  keep.beta <- matrix(0,n.samples,nc)
  keep.p    <- rep(1,n.samples)
  keep.acc  <- rep(1,n.samples)
  for(i in 1:n.samples){
    
    ## Update z 
    z <- rtruncnorm(n, a = low, b = up, mean = xbeta, sd = sqrt(1/w))
    
    ## Update w
    for (j in 1:n){
      w[j] <- rpg(2*p,z[j]-xbeta[j])
    }
    
    ## Update beta
    Omega = diag(w)
    b = prior.mn
    B = diag(prior.sd, nc)
    mulcov = solve(t(X)%*%Omega%*%X + solve(B))
    mulmean = mulcov %*% (t(X)%*%Omega%*%z+solve(B)%*%b)
    beta = mvrnorm(mu = mulmean, Sigma = mulcov)
    xbeta <- X %*% beta
    
    ## Update p through Slice sampler
    
    target = function(p){log_condi_p_y(p, y, xbeta, hyper.a , hyper.b)}
    points = sliceSample(n = n.slice, target, x.interval=c(0, 5), root.accuracy = 0.01)
    p = points[n.slice]
    
    #save the results for p and beta
    keep.p[i]    <- p
    keep.beta[i,]<- beta
    if(i%%100==0){
      print(i)
    }
  }
  # Return the posterior sample beta and p 
  list(beta=keep.beta, p=keep.p)
}

##############################################################################################
##                              MAIN PART TO BE RUN
# Generate the different simulation experiments
##############################################################################################
for(kk in 1:100){
  n = 100
  p = 3
  N = 10
  
  # Create the simulated dataset
  # 3 - dimensional
  # x         <- cbind(1,matrix(rnorm(n*2), n, 2))
  # true.beta <- c(0.5,-2,1)   #initial value of beta
  
  # 5 - dimensional
  # x         <- cbind(1,matrix(rnorm(n*4), n, 4))
  # true.beta <- c(1,-1,-3,1,3)   #initial value of beta
  
  # 10 - dimensional
  x         <- cbind(1,matrix(rnorm(n*9), n, 9))
  true.beta <- c(2.3,1,-2,1.5,-2.7,0.2,-1.4,3,-0.6,-1.2)   #initial value of beta
  
  true.p    <- exp(x%*%true.beta)/(1+exp(x%*%true.beta))
  N         <- length(true.beta)
  xbetatrue <- x%*%true.beta
  x1        <- rbeta(n,p,p)    #random beta(p,p)
  z12       <- log(x1/(1-x1))+ xbetatrue
  y         <- ifelse(z12<0,0,1)  #our new variable if z<0 then y=0 else y=1
  
  # Fit the model
  n.samples <- 5000 #number of iteration of Gibbs Sampler
  burn      <- 1000 #burn-in iteration
  n.slice   <- 4    #number of iteration within Gibbs
  prior.mn  <- 0     #prior mean for beta
  prior.sd  <- 10     #prior std deviation for beta 
  hyper.a   <- 1     #prior hyperparameters for p
  hyper.b   <- 1     #prior hyperparameters for p
  
  # Gibbs sampler function and results
  fit <- Bayes.logistic(y, x, p,true.beta, prior.mn=prior.mn, prior.sd=prior.sd,
                        hyper.a=hyper.a, hyper.b=hyper.b,n.samples=n.samples,
                        burn=burn, n.slice=n.slice)
  
  # Print the posterior mean for beta and p after a carefully burn-in:
  print(c(round(apply(fit$beta[burn:n.samples,],2,mean),3),mean(fit$p[burn:n.samples])))
}
