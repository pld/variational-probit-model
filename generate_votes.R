generate.votes <- function(true.beta.mat, votes.per.session.vec, VERBOSE=FALSE) {
  # generates the votes under the assumed voting model
  # true.beta.mat is a matrix that is num.sessions x num.object
  # each entry in that matrix is how much the person in session i likes object j
 
  votes.labels <- c("Session.ID", "Winner.ID", "Loser.ID", "Left.Choice.ID", "Right.Choice.ID")
  num.sessions <- nrow(true.beta.mat)
  num.objects <- ncol(true.beta.mat)
  num.votes <- sum(votes.per.session.vec)
  object.numbers.vec <- 1:num.objects
 
  votes <- array(NA, dim=c(num.votes, length(votes.labels)), dimnames=list(d1=NULL, d2=votes.labels))
  overall.vote.counter <- 0
  for (session in (1:num.sessions)) {
    for (vote in (1:votes.per.session.vec[session])) {
      overall.vote.counter <- overall.vote.counter + 1
      if (VERBOSE) print(paste("overall.vote.counter:", overall.vote.counter))
     
      votes[overall.vote.counter, "Session.ID"] <- session
      pair <- sample(object.numbers.vec, size=2, replace=FALSE)
      if (VERBOSE) print(paste("pair:", pair[1], " (", true.beta.mat[session, pair[1]], "),", pair[2] ," (", true.beta.mat[session, pair[2]], ")", sep=""))
      votes[overall.vote.counter, c("Left.Choice.ID", "Right.Choice.ID")] <- pair
      z <- (true.beta.mat[session, pair[1]] - true.beta.mat[session, pair[2]]) + rnorm(n=1, mean=0, sd=1)
      if (z > 0) {
        votes[overall.vote.counter, c("Winner.ID", "Loser.ID")] <- c(pair[1], pair[2])
      } else {
        votes[overall.vote.counter, c("Winner.ID", "Loser.ID")] <- c(pair[2], pair[1])
      }
      if (VERBOSE) print(votes[overall.vote.counter, ])
    }
  }
  if (any(is.na(votes))) { print(votes); stop("votes has NA"); } #check
  return(votes)
}
