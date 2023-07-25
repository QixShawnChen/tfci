

qreach <- function(x,amat,verbose = FALSE)
{
  ## Purpose: Compute possible-d-sep(x) ("pdsep")
  ## !! The non-zero entries in amat must be symmetric !!
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## - x: node of which pdsep berechnet werden soll
  ## - amat: adjacency matrix
  ##         amat[i,j] = 0 iff no edge btw i,j
  ##         amat[i,j] = 1 iff i *-o j
  ##         amat[i,j] = 2 iff i *-> j
  ## - verbose: Show checked node sequence
  ## ----------------------------------------------------------------------
  ## Author: Markus Kalisch, Date: 29 Oct 2009, 11:54
  ## Stopping:
  ## =========
  ## At every iteration, Q get's re44duced by one. It is only increased by
  ## at least one, if there are edges in A. and then, at least one
  ## edge in A. is removed. Edges are never inserted into A..
  ## Thus, either Q or A. becomes empty and the loop stops.
  ## Runtime:
  ## ========
  ## At least O(|V|), since look up in adjacency matrix is made. Assume O(|E|)>O
  ## Every edge can be visited at most twice. At each visit, there are no
  ## more than max(deg(V_i)) neighboring edges to look at. Thus, the runtime is
  ## O(2*|E| * max(deg(V_i))) = O(|E|^2) [worst case]; O(|E|) if sparse in the
  ## sense that max(deg(V_i)) is constant.
  ## Correctness:
  ## ============
  ## (A) All nodes in PSEP have a path of legal triples from x.
  ## (B) If there is a legal path from x to y in amat, at least one of them
  ## is found and y is recorded in PSEP:
  ## Suppose there is a node y != x that has a legal path from x, but is not in
  ## PSEP. y cannot be in nbrs(x), because they are added and nothing is
  ## deleted from PSEP. Hence, there must be a node z which
  ## is in PSEP but has a neighbor w in amat that has a legal path from
  ## x but is not in PSEP.
  ## Assuming that the function legal is correct, and noting that at (*) all
  ## neighbors of z in A. (tmp!) are analyzed, it follows that w is not
  ## in adj(z) in A. (but in amat). Thus, w must have been removed
  ## before from adj(z). Because of (+), it could only have been removed if
  ## u-z-w was found legal at some point. But then, w would have been added
  ## to PSEP. This is a contradiction to the assumption that w is not in PSEP.
  
  ## check: quadratic; x in V; edgemarks ok; non-zeroes symmetric
  stopifnot((ncol(amat) == nrow(amat)),x <= ncol(amat),all(amat %in% c(0,1,2)), all((amat != 0) == (t(amat != 0))))
  
  A. <- (amat != 0) ## A.[i,j] is true  <===>  edge i--j
  PSEP <- Q <- nb <- which(A.[x,])
  P <- rep.int(x, length(Q))
  A.[x,nb] <- FALSE ## delete edge to nbrs
  
  while(length(Q) > 0) {
    ## Invariants:
    ## ===========
    ## (A1) length(Q) == length(P) > 0
    ## (A2) non-zero in A. -> non-zero in amat [no non-zero added]
    ## (A3) Q[i] and P[i] are adjacent in amat [using (A2)]
    if (verbose) {
      cat("\n-------------","\n")
      cat("Queue Q:",Q,"\n")
      cat("Queue P:",P,"\n")
    }
    a <- Q[1]
    Q <- Q[-1]
    pred <- P[1] ## not empty because of (A1)
    P <- P[-1]
    if (verbose) cat("Select",pred,"towards",a,"\n")
    nb <- which(A.[a,]) ## (*)
    if (verbose) cat("Check nbrs",nb,"\nLegal:")
    
    for (b in nb) {
      ## Guaranteed: pred-a-b are a path because of (A3)
      if (lres <- legal.path(pred,a,b,amat)) {
        A.[a,b] <- FALSE ## remove b out of adj(a) in A. (+)
        Q <- c(Q,b)
        P <- c(P,a)
        PSEP <- c(PSEP,b)
      }
      if (verbose) cat(if(lres)"T" else ".", "")
    }
  }
  sort(setdiff(PSEP,x))
} ## {qreach}