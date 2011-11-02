#
# based on J. Grimmer 2010
#

# load Matrix library
library(Matrix)
library(msm)

compute.lower.bound <- function(betas, X, X.transpose.X, Y.stars, X.cov, Y, beta.mat.prior, verbose=FALSE) {
  # Compute the lower bound, expectation of logp minus expectation of logq
  #
  # Args:
  #   betas: vectors of preferences.
  #   X: Matrix of vote data, 1 for item appearing on left, -1 for item on
  #      right, otherwise 0.
  #   X.transpose.X: Cached value of t(X) %*% X.
  #   Y.stars: Current variational estimate.
  #   X.cov: Vote data covariance
  #   Y: Matrix of win data, 1 if item on left won, otherwise 0.
  #   beta.mat.prior: Covariance for the beta parameters.
  #   verbose: If TRUE, prints notes during function. Default is FALSE.
  #
  # Returns:
  #   bound: Lower bound.
  #   parts: Vector of lower bound parts.
  if (verbose) { print("calculate X^TX(\beta\beta^T + X.cov)") }
  ab <- X.transpose.X %*% (betas %*% t(betas) + X.cov)
  if (verbose) { print("calculate Tr(*)") }
  part1 <- sum(diag(ab))/2
  if (verbose) { print("calculate \beta^T\beta + Tr(X.cov)") }
  part2 <- t(betas) %*% betas + sum(diag(X.cov))
  part2 <- part2 * (beta.mat.prior[1,1])
  part2 <- as.numeric(part2/2 + (1/2)*log(det(solve(beta.mat.prior))) + (length(betas)/2*log(2*pi)))
  if (verbose) { print("calculate Y.stars^TY.stars") }
  part3 <- t(Y.stars) %*% Y.stars
  X.cov.det <- det(X.cov)
  if (verbose) { print(paste("det(X.cov): ", X.cov.det)) }
  part4 <- length(betas)/2 + (1/2)*log(X.cov.det) + (length(betas)/2)*log(2*pi)
  bound <- part1 + part2 + part3/2 + part4
  parts <- c(-part1, -part2, part3/2, part4)

  if (verbose) { print(paste("bound: ", bound)) }
  if (verbose) { print(paste("parts: ", parts)) }
  return(list(bound=bound, parts=parts))
}

variational.inference <- function(X, Y, variance=100, verbose=FALSE) {
  # Use varitational inference to estimate a preference vector.
  #
  # Args:
  #   X: Matrix of vote data, 1 for item appearing on left, -1 for item on
  #      right, otherwise 0.
  #   Y: Matrix of win data, 1 if item on left won, otherwise 0.
  #   variance: variance of beta parameters.
  #   verbose: If TRUE, prints notes during function. Default is FALSE.
  # Returns:
  #   bounds: The bounds on the last iteration after convergence.
  #   betas: The preference vector after convergence.
  #   X.cov: The vote data covariance.
  if (verbose) { print("creating beta.mat.prior") }
  # covariance for the beta parameters
  beta.mat.prior <- .symDiagonal(ncol(X), 1/variance)
  # store each updated mean vector for beta approximating distribution
  # if (verbose) { print("creating beta.var") }
  # beta.var <- matrix(NA, nrow=1000, ncol=ncol(X))

  if (verbose) { print("creating Y.stars") }
  Y.stars <- rep(0, nrow(X))
  if (verbose) { print("creating X.cov") }
  # cache X^T and X^T X
  X.transpose <- t(X)
  X.transpose.X <- X.transpose %*% X
  X.cov <- solve(X.transpose.X + beta.mat.prior)
  if (any(is.infinite(X.cov))) stop("X.cov has Inf value(s)")

  if (verbose) { print("fill Y.stars with truncated normals") }
  for (j in 1:nrow(X)) {
    Y.stars[j] <- ifelse(Y[j]==1, rtnorm(1, mean=0.5, sd=1, lower=0, upper=Inf),
      rtnorm(1, mean=-0.5, sd=1, lower=-Inf, upper=0))
  }

  # store the progress of the lower bound on the model
  bounds <- c()
  parts <- matrix(NA, nrow=1000, ncol=4)
  converged <- 0
  j <- 0
  
  # while not converged loop
  while (converged == 0) {
    j <- j + 1
    if (verbose) { print(paste("iteration j: ", j)) }

    # update beta parameters
    if (verbose) { print("update beta parameters") }
    beta.var.j <- as.vector(solve(X.transpose.X + beta.mat.prior) %*% X.transpose %*% Y.stars)
    if (length(beta.var.j) != ncol(X)) { stop("length(beta.var.j) != ncol(X)") }
    # update mean of variational dist
    if (verbose) { print("update mean of variational dist") }
    mean.var <- as.vector(X %*% beta.var.j)
    if (length(mean.var) != nrow(X)) { stop("length(mean.var) != nrow(X)") }
    mean.var.neg <- -mean.var
    
    # compute Y.stars update
    if (verbose) { print("denom E[Y.stars] update") }
    denom <- pnorm(mean.var.neg)
    if (length(denom) != nrow(X)) { stop("length(denom) != nrow(X)") }
    if (verbose) { print("numer E[Y.stars] update") }
    numer <- dnorm(mean.var.neg)
    if (length(numer) != nrow(X)) { stop("length(numer) != nrow(X)") }

    # compute E[Y.stars]
    if (verbose) { print("compute E[Y.stars] update") }
    Y.stars[which(Y==0)] <- mean.var[Y==0] + -numer[which(Y==0)]/denom[which(Y==0)]
    Y.stars[which(Y==1)] <- mean.var[Y==1] + numer[which(Y==1)]/(1 - denom[which(Y==1)])
    if (length(Y.stars) != nrow(X)) { stop("length(Y.stars) != nrow(X)") }

    # calculating lower bound
    if (verbose) { print("call lower.bound") }
    if (any(is.infinite(X.cov))) stop("X.cov has Inf value(s)")
    bounds.parts <- compute.lower.bound(beta.var.j, X, X.transpose.X, Y.stars, X.cov, Y, beta.mat.prior, verbose)
    if (verbose) { print(paste("bounds.parts$bounds: ", bounds.parts$bound)) }
    bounds[j] <- bounds.parts$bound
    if (verbose) { print(paste("bounds.parts$parts: ", bounds.parts$parts)) }
    parts[j,] <- bounds.parts$parts

    if (j > 1) {
      # check convergences
      change.in.bound <- abs(bounds[j] - bounds[j - 1])
      if (verbose) { print(paste("change.in.bound: ", change.in.bound)) }
      if (is.nan(change.in.bound)) { print("WARNING: change.in.bound is nan") }
      if ((is.nan(change.in.bound) && j > 99) || change.in.bound < 1e-8) {
        converged <- 1
        if (verbose) { print(paste("converged on iteration: ", j)) }
      }
    }
  }

  return(list(bounds=bounds, betas=beta.var.j, X.cov=X.cov))
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
    
    if (VERBOSE) { print("calling variation inference") }
    example.run <- variational.inference(X, Y, verbose=VERBOSE)
    
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
  }
  plot(mean.squared.errors, type="o")
} else {
  # load vote data
  load('fake_votes_small.RData')

  # take subset of vote data
  #votes <- votes[1:4000,]
  votes.sparse <- convert.votes.to.sparse.design.matrix(votes, VERBOSE)
  
  X <- votes.sparse$x.Matrix
  Y <- votes.sparse$y.vec
  if (VERBOSE) { print("calling variational.inference") }
  example.run <- variational.inference(X, Y, verbose=VERBOSE)
}


