#
# based on J. Grimmer 2010
#

# load Matrix library
library(Matrix)
library(msm)

# lower bound
lower.bound <- function(betas, X, Y.stars, X.cov, Y, beta.mat.prior, beta.mat.prior.half.log.det.solve) {
  ab <- t(X) %*% X %*% (betas %*% t(betas) + X.cov)
  part1 <- sum(diag(ab))/2
  part2 <- t(betas) %*% betas + sum(diag(X.cov))
  part2 <- part2 * (beta.mat.prior[1,1])
  part2 <- as.numeric(part2/2 + beta.mat.prior.half.log.det.solve + (length(betas)/2*log(2*pi)))
  part3 <- t(Y.stars) %*% Y.stars
  part4 <- length(betas)/2 + (1/2)*log(det(X.cov)) + (length(betas)/2)*log(2*pi)
  bound <- as.numeric(part1 + part2 + part3/2 + part4)
  parts <- c(-part1, -part2, part3/2, part4)

  # return bound and parts
#   bounds.parts <- list(bounds, parts)
#   names(bounds.parts) <- c('bounds', 'parts')
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
  beta.mat.prior <- .symDiagonal(ncol(X), 1/100)
  # pre-calculate as beta.mat.prior is constant
  beta.mat.prior.half.log.det.solve <- (1/2)*log(det(solve(beta.mat.prior)))

  # store each updated mean vector for beta approximating distribution
  if (VERBOSE) { print("creating beta.var") }
  beta.var <- matrix(NA, nrow=1000, ncol=ncol(X))

  if (VERBOSE) { print("creating Y.stars") }
  Y.stars <- rep(0, nrow(X))
  if (VERBOSE) { print("creating X.cov") }
  X.cov <- solve(t(X) %*% X + beta.mat.prior)

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
    beta.var.j <- solve(t(X) %*% X + beta.mat.prior) %*% t(X) %*% Y.stars

    # update mean of variational dist
    if (VERBOSE) { print("update mean of variational dist") }
    mean.var <- X %*% beta.var.j
    mean.var.neg.matrix <- as.matrix(-mean.var)
    
    # compute Y.stars update
    if (VERBOSE) { print("denom E[Y.stars] update") }
    denom <- pnorm(mean.var.neg.matrix)
    if (VERBOSE) { print("numer E[Y.stars] update") }
    numer <- dnorm(mean.var.neg.matrix)

    # compute E[Y.stars]
    if (VERBOSE) { print("compute E[Y.stars] update") }
    Y.stars[which(Y==0)] <- mean.var[Y==0] + -numer[which(Y==0)]/denom[which(Y==0)]
    Y.stars[which(Y==1)] <- mean.var[Y==1] + numer[which(Y==1)]/(1 - denom[which(Y==1)])

    # calculating lower bound
    if (VERBOSE) { print("call lower.bound") }
    bounds.parts <- lower.bound(beta.var.j, X, Y.stars, X.cov, Y, beta.mat.prior, beta.mat.prior.half.log.det.solve)
    if (VERBOSE) { print(paste("bounds.parts$bounds: ", bounds.parts$bound)) }
    bounds[j] <- bounds.parts$bound
    if (VERBOSE) { print(paste("bounds.parts$parts: ", bounds.parts$parts)) }
    parts[j,] <- bounds.parts$parts

    if (j > 1) {
      # check convergences
      change.in.bound <- abs(bounds[j] - bounds[j - 1])
      if (VERBOSE) { print(paste("change.in.bound: ", change.in.bound)) }
      if (is.nan(change.in.bound) || change.in.bound < 1e-8) {
        converged <- 1
      }
    }
  }

  # format data before returing
  data <- list(bounds, beta.var.j, X.cov)
  names(data) <- c('bounds', 'betas', 'X.cov')
  return(data)
}

# load vote data
load('fake_votes.RData')

# take subset of vote data
votes <- votes[1:4000,]

# load sparse data formatter
source('convert_votes_to_sparse_design_matrix.R')
VERBOSE <- TRUE
votes.sparse <- convert.votes.to.sparse.design.matrix(votes, VERBOSE)

X <- votes.sparse$x.Matrix
Y <- votes.sparse$y.vec

# call var inference
if (VERBOSE) { print("calling func.reg") }
example.run <- func.reg(X, Y, VERBOSE)
