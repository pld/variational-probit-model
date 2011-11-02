#
# based on J. Grimmer 2010
#

# load Matrix library
library(Matrix)
library(msm)

# lower bound
lower.bound <- function(betas, X, X.transpose.X, Y.stars, X.cov, Y, beta.mat.prior, variance) {
  if (VERBOSE) { print("calculate X^TX(\beta\beta^T + X.cov)") }
  ab <- X.transpose.X %*% (betas %*% t(betas) + X.cov)
  if (VERBOSE) { print("calculate Tr(*)") }
  part1 <- sum(diag(ab))/2
  if (VERBOSE) { print("calculate \beta^T\beta + Tr(X.cov)") }
  part2 <- t(betas) %*% betas + sum(diag(X.cov))
  part2 <- part2 * (beta.mat.prior[1,1])
  part2 <- as.numeric(part2/2 + log(det(solve(beta.mat.prior))) + (length(betas)/2*log(2*pi)))
  if (VERBOSE) { print("calculate Y.stars^TY.stars") }
  part3 <- t(Y.stars) %*% Y.stars
  # part4 is +Inf, because det(X.cov)==+Inf
  X.cov.det <- det(X.cov)
  if (VERBOSE) { print(paste("det(X.cov): ", X.cov.det)) }
  part4 <- length(betas)/2 + (1/2)*log(X.cov.det) + (length(betas)/2)*log(2*pi)
  bound <- part1 + part2 + part3/2 + part4
  parts <- c(-part1, -part2, part3/2, part4)

  # return bound and parts
  if (VERBOSE) { print(paste("bound: ", bound)) }
  if (VERBOSE) { print(paste("parts: ", parts)) }
  return(list(bound=bound, parts=parts))
}

# X is the vector of covariates
# Y is the dependent variable
func.reg <- function(X, Y, VERBOSE=FALSE) {
  # covariance for the beta parameters
  if (VERBOSE) { print("creating beta.mat.prior") }
  # create as sparse matrix
  variance <- 100
  beta.mat.prior <- .symDiagonal(ncol(X), 1/variance)
  # store each updated mean vector for beta approximating distribution
  # if (VERBOSE) { print("creating beta.var") }
  # beta.var <- matrix(NA, nrow=1000, ncol=ncol(X))

  if (VERBOSE) { print("creating Y.stars") }
  Y.stars <- rep(0, nrow(X))
  if (VERBOSE) { print("creating X.cov") }
  # cache t(X)===X^T
  X.transpose <- t(X)
  X.transpose.X <- X.transpose %*% X
  X.cov <- solve(X.transpose.X + beta.mat.prior)
  if (any(is.infinite(X.cov))) stop("X.cov has Inf value(s)")

  if (VERBOSE) { print("fill Y.stars with truncated normals") }
  for (j in 1:nrow(X)) {
    Y.stars[j] <- ifelse(Y[j]==1, rtnorm(1, mean=0.5, sd=1, lower=0, upper=Inf), rtnorm(1, mean=-0.5, sd=1, lower=-Inf, upper=0))
  }

  # store the progress of the lower bound on the model
  bounds <- c()
  parts <- matrix(NA, nrow=1000, ncol=4)
  converged <- 0
  j <- 0
  
  # while not converged loop
  while (converged == 0) {
    j <- j + 1
    if (VERBOSE) { print(paste("iteration j: ", j)) }

    # update beta parameters
    if (VERBOSE) { print("update beta parameters") }
    beta.var.j <- as.vector(solve(X.transpose.X + beta.mat.prior) %*% X.transpose %*% Y.stars)
    if (length(beta.var.j) != ncol(X)) { stop("length(beta.var.j) != ncol(X)") }
    # update mean of variational dist
    if (VERBOSE) { print("update mean of variational dist") }
    mean.var <- as.vector(X %*% beta.var.j)
    if (length(mean.var) != nrow(X)) { stop("length(mean.var) != nrow(X)") }
    mean.var.neg.matrix <- -mean.var
    
    # compute Y.stars update
    if (VERBOSE) { print("denom E[Y.stars] update") }
    denom <- pnorm(mean.var.neg.matrix)
    if (length(denom) != nrow(X)) { stop("length(denom) != nrow(X)") }
    if (VERBOSE) { print("numer E[Y.stars] update") }
    numer <- dnorm(mean.var.neg.matrix)
    if (length(numer) != nrow(X)) { stop("length(numer) != nrow(X)") }

    # compute E[Y.stars]
    if (VERBOSE) { print("compute E[Y.stars] update") }
    Y.stars[which(Y==0)] <- mean.var[Y==0] + -numer[which(Y==0)]/denom[which(Y==0)]
    Y.stars[which(Y==1)] <- mean.var[Y==1] + numer[which(Y==1)]/(1 - denom[which(Y==1)])
    if (length(Y.stars) != nrow(X)) { stop("length(Y.stars) != nrow(X)") }

    # calculating lower bound
    if (VERBOSE) { print("call lower.bound") }
    if (any(is.infinite(X.cov))) stop("X.cov has Inf value(s)")
    bounds.parts <- lower.bound(beta.var.j, X, X.transpose.X, Y.stars, X.cov, Y, beta.mat.prior, variance)
    if (VERBOSE) { print(paste("bounds.parts$bounds: ", bounds.parts$bound)) }
    bounds[j] <- bounds.parts$bound
    if (VERBOSE) { print(paste("bounds.parts$parts: ", bounds.parts$parts)) }
    parts[j,] <- bounds.parts$parts

    if (j > 1) {
      # check convergences
      change.in.bound <- abs(bounds[j] - bounds[j - 1])
      if (VERBOSE) { print(paste("change.in.bound: ", change.in.bound)) }
      if ((is.nan(change.in.bound) && j > 99) || change.in.bound < 1e-8) {
        converged <- 1
        if (VERBOSE) { print(paste("converged on iteration: ", j)) }
      }
    }
  }

  # format data before returing
  data <- list(bounds, beta.var.j, X.cov)
  names(data) <- c('bounds', 'betas', 'X.cov')
  return(data)
}

# for repeatability set rand num gen
set.seed(1)

# set verbose output
VERBOSE <- FALSE
# generate or load votes
GENERATE <- TRUE

# load sparse data formatter
source('convert_votes_to_sparse_design_matrix.R')

if (GENERATE) {
  # load vote generation code
  source('generate_votes.R')
  votes.per.session.vec <- c(1000)
  num.objects <- 9
  true.beta.mat <- t(as.matrix(runif(num.objects, -1, 1)))
  
  # run many trials 
  num.iterations <- 100
  all.betas.ordered <- matrix(NA, nrow=num.iterations, ncol=num.objects)
  mean.squared.errors <- c()
  for (iteration in 1:num.iterations) {
    print(paste("begin iteration ", iteration))
    # generate votes
    if (VERBOSE) { print(paste('generating ', votes.per.session.vec[1], ' votes')) }
    votes <- generate.votes(true.beta.mat, votes.per.session.vec, VERBOSE)
    votes.sparse <- convert.votes.to.sparse.design.matrix(votes, VERBOSE)
  
    X <- votes.sparse$x.Matrix
    Y <- votes.sparse$y.vec
    
    if (VERBOSE) { print("calling func.reg") }
    example.run <- func.reg(X, Y, VERBOSE)
    
    betas.ordered <- c()
    for (j in 1:length(votes.sparse$beta.parameters.labels)) {
      for (i in 1:length(votes.sparse$beta.parameters.identified.labels)) {
        if (votes.sparse$beta.parameters.labels[j] == votes.sparse$beta.parameters.identified.labels[i]) {
          betas.ordered <- c(betas.ordered, example.run$betas[i])
        }
      }
    }
    
    squared.errors <- (betas.ordered - as.vector(true.beta.mat))^2
    print(paste("sum of squared errors: ", sum(squared.errors)))
    print(votes.sparse$beta.parameters.labels)
    print("ordered betas:")
    print(betas.ordered)  
    print("truth:")
    print(as.vector(true.beta.mat))
    print("squared errors:")
    print(squared.errors)
    
    all.betas.ordered[iteration,] <- betas.ordered
    if (iteration > 1) {
      betas.means <- colMeans(all.betas.ordered[1:iteration,])
      print("mean over all iterations thus far:")
      print(betas.means)
      squared.errors <- (betas.means - as.vector(true.beta.mat))^2
      print(paste("sum of squared errors for means: ", sum(squared.errors)))
      mean.squared.errors <- c(mean.squared.errors, sum(squared.errors))
    } else {
      mean.squared.errors <- c(mean.squared.errors, sum(squared.errors))      
    }
    plot(mean.squared.errors, type="o")
  }
} else {
  # load vote data
  load('fake_votes_small.RData')

  # take subset of vote data
  #votes <- votes[1:4000,]
  votes.sparse <- convert.votes.to.sparse.design.matrix(votes, VERBOSE)
  
  X <- votes.sparse$x.Matrix
  Y <- votes.sparse$y.vec
  if (VERBOSE) { print("calling func.reg") }
  example.run <- func.reg(X, Y, VERBOSE)
}


