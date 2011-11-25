#
# Copyright 2011 Peter Lubell-Doughtie and the Trustees of Princeton University.
# Licensed under the BSD 3-Clause license, see file 
#
# based on J. Grimmer 2011, Political Analysis 19:32-47
#

# load Matrix library
library(Matrix)
# load msm with truncated normal
library(msm)
# start multicore libary
library(doMC)
registerDoMC()
# multicore calls run in foreach loops
library(foreach)

compute.lower.bound <- function(betas, X, X.transpose.X, Y.stars, X.var, Y,
  beta.mat.prior, verbose=FALSE, warn=TRUE) {
  # Compute the lower bound, expectation of logp minus expectation of logq
  #
  # Args:
  #   betas: vectors of preferences.
  #   X: Matrix of vote data, 1 for item appearing on left, -1 for item on
  #      right, otherwise 0.
  #   X.transpose.X: Cached value of t(X) %*% X.
  #   Y.stars: Current variational estimate.
  #   X.var: Vote data variance
  #   Y: Matrix of win data, 1 if item on left won, otherwise 0.
  #   beta.mat.prior: Variance for the beta parameters.
  #   verbose: If TRUE, prints notes during function. Default is FALSE.
  #   warn: Print warning output. Default is TRUE.
  #
  # Returns:
  #   bound: Lower bound.
  #   parts: Vector of lower bound parts.
  if (verbose) { print("calculate X^TX(\beta\beta^T + X.vra)") }
  ab <- X.transpose.X %*% (betas %*% t(betas) + X.var)
  if (verbose) { print("calculate Tr(*)") }
  part1 <- sum(diag(ab))/2
  if (verbose) { print("calculate \beta^T\beta + Tr(X.var)") }
  part2 <- t(betas) %*% betas + sum(diag(X.var))
  part2 <- part2 * (beta.mat.prior[1,1])
  part2 <- as.numeric(part2/2 + (1/2)*log(det(solve(beta.mat.prior))) +
    (length(betas)/2)*log(2*pi))
  if (verbose) { print("as.numeric part2:") }
  if (verbose) { print(part2) }
  if (verbose) { print("calculate Y.stars^TY.stars") }
  part3 <- t(Y.stars) %*% Y.stars
  X.var.det <- det(X.var)
  if (verbose) { print(paste("det(X.var): ", X.var.det)) }
  part4 <- length(betas)/2 + (1/2)*log(X.var.det) + (length(betas)/2)*log(2*pi)
  bound <- part1 + part2 + part3/2 + part4
  parts <- c(-part1, -part2, part3/2, part4)

  if (verbose) { print(paste("bound: ", bound)) }
  if (verbose) { print(paste("parts: ", parts)) }
  return(list(bound=bound, parts=parts))
}

variational.inference <- function(X, Y, variance=100, verbose=FALSE, warn=TRUE,
  max.iterations=1000) {
  # Use varitational inference to estimate a preference vector.
  #
  # Args:
  #   X: Matrix of vote data, 1 for item appearing on left, -1 for item on
  #      right, otherwise 0.
  #   Y: Matrix of win data, 1 if item on left won, otherwise 0.
  #   variance: Variance of beta parameters. Default is 100.
  #   verbose: If TRUE, prints notes during function. Default is FALSE.
  #   warn: Print warning output. Default is TRUE.
  #   max.iterations: Maximum iterations before assuming convergence. Default
  #      is 1000.
  #
  # Returns:
  #   bounds: The bounds on the last iteration after convergence.
  #   betas: The preference vector after convergence.
  #   X.var: The vote data variance.
  if (verbose) { print("creating beta.mat.prior") }
  # variance for the beta parameters
  # no functional difference between sparse and non-sparse diagonal matrices
  beta.mat.prior <- .symDiagonal(ncol(X), 1/variance)

  if (verbose) { print("creating X.var") }
  # cache X^T, X^T X, (X^T X + \beta)^{-1}, (X^T X + \beta)^{-1} X^T
  X.transpose <- t(X)
  X.transpose.X <- X.transpose %*% X
  X.var <- solve(X.transpose.X + beta.mat.prior)
  if (any(is.infinite(X.var))) stop("X.var has Inf value(s)")
  X.var.X.transpose <- X.var %*% X.transpose

  if (verbose) { print("creating Y.stars") }
  Y.stars <- rep(0, nrow(X))
  if (verbose) { print("fill Y.stars with truncated normals") }
  for (j in 1:nrow(X)) {
    Y.stars[j] <- ifelse(Y[j]==1, rtnorm(1, mean=0.5, sd=1, lower=0, upper=Inf),
      rtnorm(1, mean=-0.5, sd=1, lower=-Inf, upper=0))
  }

  # store each updated mean vector for beta approximating distribution
  if (verbose) { print("creating beta.var") }
  beta.var <- matrix(NA, nrow=max.iterations, ncol=ncol(X))
  # store the progress of the lower bound on the model
  bounds <- c()
  parts <- matrix(NA, nrow=max.iterations, ncol=4)
  converged <- 0
  j <- 0 
  
  # while not converged loop
  while (converged == 0 && j < max.iterations) {
    j <- j + 1
    if (verbose) { print(paste("iteration j: ", j)) }

    # update beta parameters
    if (verbose) { print("update beta parameters") }
    beta.var[j, ] <- as.vector(X.var.X.transpose %*% Y.stars)
    if (length(beta.var[j, ]) != ncol(X)) {
      stop("length(beta.var[j, ]) != ncol(X)")
    }
    # update mean of variational dist
    if (verbose) { print("update mean of variational dist") }
    mean.var <- as.vector(X %*% beta.var[j, ])
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
    Y.stars[which(Y==0)] <- mean.var[Y==0] +
      -numer[which(Y==0)]/denom[which(Y==0)]
    Y.stars[which(Y==1)] <- mean.var[Y==1] +
      numer[which(Y==1)]/(1 - denom[which(Y==1)])
    if (length(Y.stars) != nrow(X)) { stop("length(Y.stars) != nrow(X)") }

    # calculating lower bound
    if (verbose) { print("call lower.bound") }
    if (any(is.infinite(X.var))) stop("X.var has Inf value(s)")
    bounds.parts <- compute.lower.bound(beta.var[j, ], X, X.transpose.X, 
      Y.stars, X.var, Y, beta.mat.prior, verbose)
    if (verbose) { print(paste("bounds.parts$bounds: ", bounds.parts$bound)) }
    bounds[j] <- bounds.parts$bound
    if (verbose) { print(paste("bounds.parts$parts: ", bounds.parts$parts)) }
    parts[j, ] <- bounds.parts$parts

    if (j > 1) {
      # check convergences
      change.in.bound <- abs(bounds[j] - bounds[j - 1])
      if (verbose) { print(paste("change.in.bound: ", change.in.bound)) }
      if (warn && is.nan(change.in.bound)) {
        print("WARNING: change.in.bound is nan")
      }
      if (is.nan(change.in.bound) || change.in.bound < 1e-8) {
        converged <- 1
        if (verbose) { print(paste("converged on iteration: ", j)) }
      }
    }
  }
  if (warn && j >= max.iterations) {
    print(paste("WARNING: exceeded max.iterations ", max.iterations))
    print(bounds[(max.iterations - 10):max.iterations])
  }

  return(list(bounds=bounds, betas=beta.var[j, ], X.var=X.var))
}

order.betas <- function(votes.sparse, estimated.beta) {
  betas.ordered <- c()
  labels <- votes.sparse$beta.parameters.identified.labels
  for (j in 1:length(votes.sparse$beta.parameters.labels)) {
    for (i in 1:length(labels)) {
      if (votes.sparse$beta.parameters.labels[j] == labels[i]) {
        betas.ordered[j] <- estimated.beta[i]
        break
      }
    }
    if (is.null(betas.ordered[j])) {
      print("parameters.identified.labels")
      print(labels)
      print("estimated.beta")
      print(estimated.beta)
      print("parameters.labels")
      print(votes.sparse$beta.parameters.labels)
      print("betas.ordered")
      print(betas.ordered)        
      stop(paste("no ordered beta assigned to j: ", j))
    }
  }
  return(list(betas.ordered=betas.ordered))
}
test.variational.inference <- function(warn=TRUE, verbose=FALSE, generate=TRUE,
  change.betas=TRUE) {
  # Run tests of varitational inference.
  #
  # Args:
  #   warn: Print warning output. Default is TRUE.
  #   verbose: If TRUE, prints notes during function. Default is FALSE.
  #   generate: Generate votes or load. Default is TRUE.
  #   change.betas: Generate new true betas on each replication.
  #
  # Returns:
  #   nothing
  
  # load sparse data formatter
  source('convert_votes_to_sparse_design_matrix.R')
  
  if (generate) {
    # load vote generation code
    source('generate_votes.R')
    
    # run many trials 
    num.replications <- 1500
    sample.size.vec <- seq(1000, 110000, 10000)

    means.per.sizes <- matrix(nrow=length(sample.size.vec),
      ncol=num.replications)
    
    for (votes.size in 1:length(sample.size.vec)) {
      num.votes <- sample.size.vec[votes.size]
      votes.per.session.vec <- c(num.votes)
      print(paste("vote size ", votes.size, " will generate ", num.votes,
        " votes"))

      num.objects <- 9
      true.beta.mat <- t(as.matrix(runif(num.objects, -2, 2)))
      
      all.betas.ordered <- matrix(NA, nrow=num.replications, ncol=num.objects)
      
      all.mse <- c()
      all.mse <- foreach (rep = 1:num.replications, .combine = "c") %dopar% {
        print(paste("replication ", rep))
        # generate votes
        if (verbose) {
          print(paste('generating ', num.votes, ' votes'))
        }
        votes <- generate.votes(true.beta.mat, votes.per.session.vec, verbose)
        votes.sparse <- convert.votes.to.sparse.design.matrix(votes, verbose)
      
        X <- votes.sparse$x.Matrix
        Y <- votes.sparse$y.vec
        
        if (verbose) { print("calling variational inference") }
        results.var <- variational.inference(X, Y, verbose=verbose, warn=warn)
        betas <- results.var$betas
        
        if (verbose) {
          iter.mse <- c()
          for (j in 1:nrow(betas)) {
            betas.ordered <- order.betas(votes.sparse, betas[j, ])$betas.ordered
            iter.mse[j] <- log(sum(betas.ordered - 
              as.vector(true.beta.mat))^2/ncol(betas))
          }
          plot(iter.mse, main=paste("Log MSE vs Iteration Vote Size ",
            num.votes, " Replication ", rep))
        }

        estimated.beta <- betas
        betas.ordered <- order.betas(votes.sparse, estimated.beta)$betas.ordered
        
        squared.errors <- (betas.ordered - as.vector(true.beta.mat))^2
        mse <- sum(squared.errors)/length(squared.errors)
        if (verbose) {
          print(paste("mse: ", mse))
          print(votes.sparse$beta.parameters.labels)
          print("ordered betas:")
          print(betas.ordered)  
          print("truth:")
          print(as.vector(true.beta.mat))
          print("squared errors:")
          print(squared.errors)
        }
        
        all.betas.ordered[rep, ] <- betas.ordered
#         all.mse[rep] <- mse
        if (change.betas) {
          true.beta.mat <- t(as.matrix(runif(num.objects, -1, 1)))
        }
        mse
      }
      # end loop over replications
      means.per.sizes[votes.size, ] <- all.mse
      if (verbose) {
        plot(all.mse, type="o", main=paste(num.votes, " Votes"),
          xlab="Replication", ylab="Mean Sum of Squares Error")
      }
      cumulative.mean.all.mse <- c()
      for (j in 1:length(all.mse)) {
        cumulative.mean.all.mse[j] <- mean(all.mse[1:j])
      }
      # we expect the cumulative mean to stabilize over more replications
      plot(cumulative.mean.all.mse, type="o", main=paste(num.votes, " Votes"),
        xlab="Replication", ylab="Cumulative Mean Sum of Squares Error")
    }
    # end loop over vote sizes
    indices <- matrix(ncol=ncol(means.per.sizes), nrow=nrow(means.per.sizes))
    for (j in 1:nrow(means.per.sizes)) {
      indices[j, ] <- rep(j, ncol(means.per.sizes))
    }
    # we expect the MSE to decrease over larger vote sizes
    plot(indices, means.per.sizes, xlab="Number of Votes",
      ylab="Mean Sum of Squares Error", xaxt="n",
      main=paste(num.replications, " replications per size"))
    axis(1, c(1:length(sample.size.vec)), sample.size.vec)
    for (j in 1:nrow(means.per.sizes)) {
      points(j, mean(means.per.sizes[j, ]), col="red", pch=18, cex=3)
    }
    return (list(means.per.sizes=means.per.sizes))
  } else {
    # load vote data
    load('fake_votes_small.RData')
  
    # take subset of vote data
    votes.max <- 0
    if (votes.max > 0) { votes <- votes[1:votes.max, ] }
    votes.sparse <- convert.votes.to.sparse.design.matrix(votes, verbose)
    
    X <- votes.sparse$x.Matrix
    Y <- votes.sparse$y.vec
    if (verbose) { print("calling variational.inference") }
    example.run <- variational.inference(X, Y, verbose=verbose, warn=warn)
    print(example.run$betas)
    print(votes.sparse$beta.parameters.identified.labels)
  }
}

# for repeatability set rand num gen
set.seed(1)
# run tests
test.results <- test.variational.inference()
