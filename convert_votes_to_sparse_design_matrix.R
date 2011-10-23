library(Matrix)

get.objects <- function(labels, VERBOSE=FALSE) {
  # takes a bunch of labels
  # and returns the objects in those labels
  # Example input: c("o.1023.s242354", "o.123343.s1234234")
  # Example output: c("o.1023", "o.123343")
   
  pieces.list <- strsplit(labels, "\\.");
  if (length(pieces.list[[1]]) != 4) stop("ERROR: length of pieces list is not correct");
  object.numbers.temp <- sapply(pieces.list, "[[", 2);
  # for more on this command see: https://stat.ethz.ch/pipermail/r-help/2004-September/057550.html
  results <- paste("o.", object.numbers.temp, sep="");
 
  if (length(results) != length(labels)) stop("ERROR: input and output not same size");
  return(results);
}

get.object <- function(parameter.str) {
  split.parameter.str <- strsplit(parameter.str, split="\\.");
  result <- paste(split.parameter.str[[1]][1], ".", split.parameter.str[[1]][2], sep="");
  return(result);
}

convert.votes.to.sparse.design.matrix <- function(set.of.votes, VERBOSE=FALSE) {
  # takes a set of votes and turns it into a design matrix for the probit model
  if (VERBOSE==TRUE) print("Converting votes file to full design matrix for probit . . . ");
  object.numbers <- sort(unique(c(set.of.votes[,"Winner.ID"], set.of.votes[,"Loser.ID"])));
  object.labels <- paste("o.", object.numbers, sep="");
  session.numbers <- sort(unique(set.of.votes[,"Session.ID"]));
  session.labels <- paste("s.", session.numbers, sep="");
  if (VERBOSE==TRUE) { print("head(object.numbers):"); print(head(object.numbers)); }
  if (VERBOSE==TRUE) { print("head(session.numbers):"); print(head(session.numbers)); }

  # check input
  if (any(set.of.votes[, "Left.Choice.ID"]==set.of.votes[,"Right.Choice.ID"])) stop("ERROR: Some objects seem to compete against themself.  That's not possible.");
  if (!all(object.numbers %in% set.of.votes[, "Winner.ID"])) stop("ERROR: not every object has won"); 
  if (!all(object.numbers %in% set.of.votes[, "Loser.ID"])) stop("ERROR: not every object has lost");
  
  combo.df <- expand.grid(session.labels, object.labels);
  beta.parameters.labels <- paste(combo.df[,2], ".", combo.df[,1], sep="");
  # the order of beta parameters labels matches how they would fill 
  # matrix with num.session rows and num.objects columns
  # e.g., "o.1.s.1", "o.1.s.2" "o.1.s.3", "o.2.s.1", "o.2.s.2"
  if (VERBOSE==TRUE) { print("head(beta.parameters.labels):"); print(head(beta.parameters.labels)); }
  num.votes <- nrow(set.of.votes);
  num.sessions <- length(session.labels);
  
  # design.matrix, we will create it directly in triplet form
  # need, i.vec, j.vec, and x.vec
  # i.vec will be 1, 2, ... num.votes, 1, 2,  ... num.votes
  # j.vec (column of 1st object.session on left, 2nd object.session on left . . . .1st object session on right
  # x.vec will be 1, 1, ... 1, -1, -1, ... -1
  
  i.vec <- rep((1:num.votes), 2);
  if (VERBOSE) print("i.vec complete");
  
  left.object.session.vec <- paste("o.", as.character(set.of.votes[,"Left.Choice.ID"]), ".s.", as.character(set.of.votes[,"Session.ID"]), sep="");
  right.object.session.vec <- paste("o.", as.character(set.of.votes[,"Right.Choice.ID"]), ".s.", as.character(set.of.votes[,"Session.ID"]), sep="");
  j.vec <- c(match(left.object.session.vec, beta.parameters.labels), match(right.object.session.vec, beta.parameters.labels))
  if (VERBOSE) print("j.vec complete");

  x.vec <- c(rep(1, num.votes), rep(-1, num.votes));
  if (VERBOSE) print("x.vec complete");

  x.Matrix.full <- sparseMatrix(i=i.vec, j=j.vec, x=x.vec, dims=c(num.votes, length(beta.parameters.labels)), dimnames=list(d1=NULL, d2=beta.parameters.labels), index1=TRUE);

  if (VERBOSE) print("Beginning finding identified and unidentified parameters")
  beta.parameters.identified.labels <- unique(c(left.object.session.vec, right.object.session.vec));
  if (VERBOSE) { print("head(beta.parameters.identified.labels):"); print(head(beta.parameters.identified.labels)); }
  beta.parameters.not.identified.indexes <- !(beta.parameters.labels %in% beta.parameters.identified.labels);
  beta.parameters.not.identified.labels <- beta.parameters.labels[beta.parameters.not.identified.indexes];
  if (VERBOSE) { print("head(beta.parameters.not.identified.labels):"); print(head(beta.parameters.not.identified.labels)); }
  beta.parameters.not.identified.object.labels <- get.objects(beta.parameters.not.identified.labels, VERBOSE=FALSE);
  if (VERBOSE) { print("head(beta.parameters.not.identified.obejct.labels):"); print(head(beta.parameters.not.identified.object.labels)); }
  beta.parameters.identified.object.labels <- get.objects(beta.parameters.identified.labels, VERBOSE=FALSE);
  if (VERBOSE) { print("head(beta.parameters.identified.object.labels):"); print(head(beta.parameters.identified.object.labels)); }

  if (length(beta.parameters.identified.labels)
      + length(beta.parameters.not.identified.labels)
      != length(beta.parameters.labels)) stop("ERROR with length of beta.parameter vectors")
  if (length(beta.parameters.not.identified.labels)
      != length(beta.parameters.not.identified.labels)) stop("ERROR strip beta.parameters.not.identified")
  if (any(is.na(beta.parameters.identified.labels))) stop("ERROR beta.parameters.identified.labels has NA");  
  if (any(is.na(beta.parameters.not.identified.labels))) stop("ERROR beta.parameters.not.identified.labels has NA");  
  if (any(is.na(beta.parameters.not.identified.object.labels))) stop("ERROR beta.parameters.not.identified.object.labels has NA");  
  if (VERBOSE) print("Finished finding identified and unidentified parameters")

  if (VERBOSE) { 
    print(paste("beta.parameters:", length(beta.parameters.labels)));
    print(paste("beta.parameters.identified:", length(beta.parameters.identified.labels)));
    print(paste("percent of beta parameters identified:",  100*round(length(beta.parameters.identified.labels)/length(beta.parameters.labels), 2)));
  }
  if (VERBOSE) print("create x.Matrix")
  x.Matrix <- x.Matrix.full[, beta.parameters.identified.labels];
  if (ncol(x.Matrix) != length(beta.parameters.identified.labels)) stop("ERROR: problem creating x.Matrix");
  if (VERBOSE) print("create x.Matrix.aug")
  x.Matrix.aug <- rBind(x.Matrix, Diagonal(x=rep(1, ncol(x.Matrix))));

  # y.vec (1 if object on left wins, 0 if object on right wins)
  if (VERBOSE) print("create y.vec")
  y.vec <- as.numeric(set.of.votes[,"Winner.ID"]==set.of.votes[,"Left.Choice.ID"]);
    
  # check output
  if (!(all(y.vec %in% c(0,1)))) stop("ERROR: y.vec contains values other than 0 and 1");
  if (any(is.na(y.vec))) stop("ERROR: y.vec contains NA");
  if (!(all(rowSums(x.Matrix.full) == 0))) stop("ERROR: rowSums in x.Matrix.full are not 0");
  if (!(all(rowSums(x.Matrix) == 0))) stop("ERROR: rowSums in x.Matrix are not 0");
  if (!(all(x.Matrix.full@x %in% c(-1, 0, 1)))) stop("ERROR: x.Matrix.full has values other than -1, 0, 1");
  if (!(all(x.Matrix@x %in% c(-1, 0, 1)))) stop("ERROR: x.Matrix has values other than -1, 0, 1");
  if (any(is.na(x.Matrix.full))) stop("ERROR: x.Matrix.full has NA");
  if (any(is.na(x.Matrix))) stop("ERROR: x.Matrix has NA");
  return(list(y.vec=y.vec, 
              x.Matrix.full=x.Matrix.full, 
              x.Matrix=x.Matrix,
              x.Matrix.aug=x.Matrix.aug,
              beta.parameters.labels=beta.parameters.labels,
              beta.parameters.identified.labels=beta.parameters.identified.labels,
              beta.parameters.identified.object.labels=beta.parameters.identified.object.labels,
              beta.parameters.not.identified.labels=beta.parameters.not.identified.labels,
              beta.parameters.not.identified.object.labels=beta.parameters.not.identified.object.labels,
              object.labels=object.labels, object.numbers=object.numbers,
              session.labels=session.labels, session.numbers=session.numbers,
              num.sessions=num.sessions, num.votes=num.votes));
  
}

