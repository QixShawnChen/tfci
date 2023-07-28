

## only called in fci() [by default:  doPdsep=TRUE]
tpdsep <- function (skel, suffStat, indepTest, p, sepset, alpha, pMax, m.max = Inf,
                    pdsep.max = Inf, NAdelete = TRUE, unfVect = NULL,
                    #biCC = FALSE, 
                    fixedEdges = NULL, 
                    verbose = FALSE,
                    tiers = NULL) ## FIXME: verbose : 2 --> qreach(verbose)
{
  ## Purpose: Compute Possible-D-SEP for each node, perform the conditional
  ##          independent tests and adapt graph accordingly
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## - skel: Graph object returned by function skeleton   
  ## - suffStat, indepTest: info for the independence tests
  ## - p: number of nodes in the graph
  ## - sepset: Sepset that was used for finding the skeleton
  ## - alpha: niveau for the tests
  ## - pMax: Maximal p-values during estimation of skeleton
  ## - m.max: maximal size of the conditioning sets
  ## - pdsep.max: maximaum size of conditioning set for Possible-D-SEP
  ## - unfVect: vector containing the unfaithful triples, used for the
  ##   conservative orientation of the v-structures
  ## - biCC: if the biconnected components have to be used
  ## - fixedEdges: Edges marked here are not changed (logical)
  ## ----------------------------------------------------------------------
  ## Value:
  ## - G: Updated boolean adjacency matrix
  ## - sepset: Updated sepsets
  ## - pMax: Updated pMax
  ## - allPdsep: Possible d-sep for each node [list]
  ## ----------------------------------------------------------------------
  ## Author: Markus Kalisch, Date:  9 Dec 2009
  ## Modification: Diego Colombo; Martin Maechler; Joris Mooij
  
  ## cat("NOW ENTERING TPDSEP...", "\n")
  
  ##change1 starts
  ################################################# adding the step to guarantee the input of tiers parameter is valid.
  ## if no tiers are specified, everything is tier 0
  if (is.null(tiers)) {
    tiers <- rep(0, p)
  } else {
    ## check if 'tiers' are correctly specified
    if (!is.numeric(tiers)) {stop("'tiers' must be a numeric vector")}
    if (length(tiers) != p) {stop("length of 'tiers' does not match 'p' or length of 'labels'")}
  }
  #################################################
  ##change1 ends
  
  G <- (as(skel, "matrix") != 0)
  n.edgetests <- rep(0, 1000)
  ord <- 0L
  allPdsep.tmp <- vector("list", p)   
  
  #if(biCC)
  #  conn.comp <- lapply(biConnComp(skel), as.numeric)
  
  if (any(is.null(fixedEdges))) { ## MM: could be sparse
    fixedEdges <- matrix(rep(FALSE, p * p), nrow = p, ncol = p)
  }
  else if (!identical(dim(fixedEdges), c(p, p)))
    stop("Dimensions of the dataset and fixedEdges do not agree.")
  else if (!identical(fixedEdges, t(fixedEdges)) )
    stop("fixedEdges must be symmetric")
  
  if (any(G)) {
    amat <- G
    ind <- which(G, arr.ind = TRUE)       
    storage.mode(amat) <- "integer" # (TRUE, FALSE) -->  (1, 0)
    #-------------change2 starts
    ### INSERT 2 AT amat[j,i] WHEREVER THERE IS AN EDGE CONNECTING FUTURE at j TO PAST at i ###
    for(i in 1:p) {
      for(j in i:p) {
        if(amat[i,j]==1) cat("tier i = ", tiers[i]," tier j = ", tiers[j], "\n")
        if(amat[i,j] == 1 && tiers[j] > tiers[i]) { 
          amat[i,j] <- 2
          amat[j,i] <- 1
        } }
    }
    #------------change2 ends
    
    ## Orient colliders
    if (verbose) cat("\nCompute collider:\n")
    for (i in seq_len(nrow(ind))) {
      x <- ind[i, 1]
      y <- ind[i, 2]
      allZ <- setdiff(which(amat[y, ] != 0), x)
      for (z in allZ) {
        if (amat[x, z] == 0 &&
            !(y %in% sepset[[x]][[z]] ||
              y %in% sepset[[z]][[x]])) {
          
          if (length(unfVect) == 0) { ## normal version -------------------
            amat[x, y] <- amat[z, y] <- 2
            if (verbose) cat("\n",x,"*->", y, "<-*", z, "\n")
          }
          else { ## conservative version : check if x-y-z is faithful
            if (!any(unfVect == triple2numb(p,x,y,z), na.rm = TRUE) &&
                !any(unfVect == triple2numb(p,z,y,x), na.rm = TRUE)) {
              amat[x, y] <- amat[z, y] <- 2
              if (verbose)
                cat("\n",x,"*->", y, "<-*", z, "\n")
            }
          }
        }
      } ## for( z )
    } ## for( i  )
    
    allPdsep <- lapply(1:p, qreach, amat = amat)# verbose = (verbose >= 2) ## NEED TO ORIENT FUTURE/PAST + COLLIDERS BEFORE ENTERING THIS STEP ##
    allPdsep.tmp <- vector("list", p)
    for(x in seq_len(p)) {
      if(verbose) cat("\nPossible D-Sep of", x, "is:", allPdsep[[x]], "\n")
      if (any(an0 <- amat[x, ] != 0)) {
        tf1 <- setdiff(allPdsep[[x]], x)
        adj.x <- which(an0)
        for (y in adj.x)
          if( !fixedEdges[x,y] ) {
            if(verbose) cat(sprintf("\ny = %3d\n.........\n", y))
            tf <- setdiff(tf1, y)
            diff.set <- setdiff(tf, adj.x)
            
            
            ##------------ change3 starts
            
            tiers.diff.set = tiers
            vertices_vec = seq(1, p)
            vertices_non = setdiff(vertices_vec, diff.set)  ##the step is to find vertices not in diff.set
            for(i in vertices_non) {  ##the step is to change the 'tiers' value(s) of the non-existing element to -1
              tiers.diff.set[i] = -Inf
            }
            
            if(all(sapply(tiers.diff.set, is.nan))) {
              next       ##if all elements are -1, that means all the vertices in conditioning set are in the future.
            }
            #need to consider what if all the values in vector tiers.diff.set are NULL.
            
            ##------------ change3 ends
            
            
            ## bi-connected components
            #            if (biCC) {
            #              for(cci in conn.comp) {
            #                if (x %in% cci && y %in% cci)
            #                  break ## found it
            #            }
            #             bi.conn.comp <- setdiff(cci, c(x,y))
            #              tf <- intersect(tf, bi.conn.comp)
            #              if (verbose) {
            #                cat("There is an edge between",x,"and",y,"\n")
            #                cat("Possible D-Sep of", x,
            #                    "intersected with the biconnected component of",x,"and",y,
            #                   "is:", tf, "\n")
            #            }
            #          } ## if(biCC)
            
            
            
            allPdsep.tmp[[x]] <- c(tf,y) ## you must add y to the set
            ## for the large scale simulations, we need to stop the algorithm if
            ## it takes to much time, i.e. sepset>25
            if (length(tf) > pdsep.max) {
              if(verbose)
                cat("Size of Possible-D-SEP bigger than",pdsep.max,
                    ". Break the search for the edge between", x,"and",y,"\n")
            } else if (length(diff.set) > 0) {
              done <- FALSE
              ord <- 0L
              while (!done && ord < min(length(tf), m.max)) {
                ord <- ord + 1L
                if(verbose) cat("ord = ", ord, "\n")
                if (ord == 1) {
                  for (S in diff.set) {
                    
                    ##---------------------  change5—1 starts
                    
                    if(max(tiers.diff.set[[S]]) > max(tiers[x], tiers[y])) { #skip the situation when max tiers in diff.set is greater than the maximum of tiers values of vertex X and vertex Y.
                      next
                    }
                    
                    ##---------------------  change5-1 ends
                    
                    pval <- indepTest(x, y, S, suffStat)
                    n.edgetests[ord + 1] <- n.edgetests[ord + 1] + 1
                    if (is.na(pval))
                      pval <- as.numeric(NAdelete) ## = if(NAdelete) 1 else 0
                    if (pval > pMax[x, y])
                      pMax[x, y] <- pval
                    if (pval >= alpha) {
                      amat[x, y] <- amat[y, x] <- 0
                      sepset[[x]][[y]] <- sepset[[y]][[x]] <- S
                      done <- TRUE
                      if (verbose)
                        cat("x=", x, " y=", y, " S=", S, ": pval =", pval, "\n")
                      break
                    }
                  }
                }
                else { ## ord > 1
                  tmp.combn <- combn(tf, ord) ## has  choose( |tf|, ord ) columns
                  if (ord <= length(adj.x)) {
                    for (k in seq_len(ncol(tmp.combn))) {
                      S <- tmp.combn[, k]
                      if (!all(S %in% adj.x)) {
                        n.edgetests[ord + 1] <- n.edgetests[ord + 1] + 1
                        
                        ##---------------------  change5-2 starts
                        
                        if(max(tiers.diff.set[[S]]) > max(tiers[x], tiers[y])) { #skip the situation when max tiers in diff.set is                               greater than the maximum of tiers values of vertex X and vertex Y.
                          next
                        }
                        
                        ##---------------------  change5-2 ends
                        
                        pval <- indepTest(x, y, S, suffStat)
                        if (is.na(pval))
                          pval <- as.numeric(NAdelete) ## = if(NAdelete) 1 else 0
                        if(pMax[x, y] < pval)
                          pMax[x, y] <- pval
                        if (pval >= alpha) {
                          amat[x, y] <- amat[y, x] <- 0
                          sepset[[x]][[y]] <- sepset[[y]][[x]] <- S
                          done <- TRUE
                          if (verbose)
                            cat("x=", x, " y=", y, " S=", S, ": pval =", pval, "\n")
                          break
                        }
                      }
                    } ## for(k ..)
                  }
                  else { ## ord > |adj.x| :
                    ## check all combinations; no combination has been tested before
                    for (k in seq_len(ncol(tmp.combn))) {
                      S <- tmp.combn[, k]
                      n.edgetests[ord + 1] <- n.edgetests[ord + 1] + 1
                      
                      ##---------------------  change5-3 starts
                      
                      if(max(tiers.diff.set[[S]]) > max(tiers[x], tiers[y])) { #skip the situation when max tiers in diff.set is                           greater than the maximum of tiers values of vertex X and vertex Y.
                        next
                      }
                      
                      ##---------------------  change5-3 ends
                      
                      pval <- indepTest(x, y, S, suffStat)
                      if (is.na(pval))
                        pval <- as.numeric(NAdelete) ## = if(NAdelete) 1 else 0
                      if(pMax[x, y] < pval)
                        pMax[x, y] <- pval
                      if (pval >= alpha) {
                        amat[x, y] <- amat[y, x] <- 0
                        sepset[[x]][[y]] <- sepset[[y]][[x]] <- S
                        done <- TRUE
                        if (verbose)
                          cat("x=", x, " y=", y, " S=", S, ": pval =", pval, "\n")
                        break
                      }
                    } ## for(k ..)
                  } ## else: { ord > |adj.x| }
                } ## else
                
              } ## while(!done ..)
            }
            
          } ## for(y ..)
        
      } ## if(any( . ))
      
    } ## for(x ..)
    G[amat == 0] <- FALSE
    G[amat == 1] <- TRUE
    G[amat == 2] <- TRUE
    
  } ## if(any(G))
  
  list(G = G, sepset = sepset, pMax = pMax, allPdsep = allPdsep.tmp,
       max.ord = ord, n.edgetests = n.edgetests[1:(ord + 1)])
} ## {tpdsep}
