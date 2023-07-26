library(pcalg)

source("tcheckTriple.R")
source("tpc.cons.intern.R")
source("tpdsep.R")
source("tskeleton.R")

tfci <- function(suffStat, indepTest, alpha, labels, p,
                 skel.method = c("stable", "original"),
                 type = c("normal", "anytime", "adaptive"),
                 fixedGaps = NULL, fixedEdges = NULL, NAdelete = TRUE,
                 m.max = Inf, pdsep.max = Inf, rules = rep(TRUE, 10),
                 doPdsep = TRUE, 
                 biCC = FALSE, conservative = FALSE, maj.rule = TRUE, 
                 selectionBias = TRUE,
                 jci = c("0","1","12","123"), contextVars = NULL, 
                 verbose = FALSE, 
                 tiers = NULL, context.all = NULL, context.tier = NULL)
{
  ## For FCI:
  ## Purpose: Perform FCI-Algorithm, i.e., estimate PAG
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## - suffStat, indepTest: info for the tests
  ## - p: number of nodes in the graph
  ## - alpha: Significance level of individual partial correlation tests
  ## - verbose: 0 - no output, 1 - detailed output
  ## - fixedGaps: the adjacency matrix of the graph from which the algorithm
  ##      should start (logical); gaps fixed here are not changed
  ## - fixedEdges: Edges marked here are not changed (logical)
  ## - NAdelete: delete edge if pval=NA (for discrete data)
  ## - m.max: maximum size of conditioning set
  ## - pdsep.max: maximaum size of conditioning set for Possible-D-SEP
  ## - rules: array of length 10 wich contains TRUE or FALSE corresponding
  ##          to each rule. TRUE means the rule will be applied.
  ##          If rules==FALSE only R0 (minimal pattern) will be used
  ## - doPdsep: compute possible dsep
  ## - biCC: TRUE or FALSE variable containing if biconnected components are
  ##         used to compute pdsep
  ## - conservative: TRUE or FALSE defining if
  ##          the v-structures after the pdsep
  ##          have to be oriented conservatively or nor
  ## - maj.rule: TRUE or FALSE variable containing if the majority rule is
  ##             used instead of the normal conservative
  ## - labels: names of the variables or nodes
  ## - type: it specifies the version of the FCI that has to be used.
  ##         Per default it is normal, the normal FCI algorithm. It can also be
  ##         anytime for the Anytime FCI and in this cas m.max must be specified;
  ##         or it can be adaptive for Adaptive Anytime FCI and in this case
  ##         m.max must not be specified.
  ## - numCores: handed to skeleton(), used for parallelization
  ## - selectionBias: allow for selection bias (default: TRUE)
  ## - jci: specifies the JCI background knowledge that is used; can be either:
  ##     "0"   no JCI background knowledge (default),
  ##     "1"   JCI assumption 1 only (i.e., no system variable causes any context variable),
  ##     "12"  JCI assumptions 1 and 2 (i.e., no system variable causes any context variable,
  ##           and no system variable is confounded with any context variable),
  ##     "123" all JCI assumptions 1, 2 and 3
  ## - contextVars: subset of variable indices that will be treated as context variables
  ## ----------------------------------------------------------------------
  ## Author: Markus Kalisch, Date: Dec 2009; update: Diego Colombo, 2012; Martin Maechler, 2013; Joris Mooij, 2020
  
  cl <- match.call()
  if(!missing(p)) stopifnot(is.numeric(p), length(p <- as.integer(p)) == 1, p >= 2)
  if(missing(labels)) {
    if(missing(p)) stop("need to specify 'labels' or 'p'")
    labels <- as.character(seq_len(p))
  } else { ## use labels ==> p  from it
    stopifnot(is.character(labels))
    if(missing(p)) {
      p <- length(labels)
    } else if(p != length(labels))
      stop("'p' is not needed when 'labels' is specified, and must match length(labels)")
    else
      message("No need to specify 'p', when 'labels' is given")
  }
  
  if (is.null(tiers)) {
    ## if no tiers are specified, everything is tier 1
    tiers <- rep(1, p)
  } else {
    ## check if 'tiers' are correctly specified
    if (!is.numeric(tiers)) {stop("'tiers' must be a numeric vector")}
    if (length(tiers) != p) {stop("length of 'tiers' does not match 'p' or length of 'labels'")}
  }
  
  
  
  
  ##Change1 starts
  
  if (!is.null(context.all)) {
    if (is.character(context.all)) {
      if (!all(context.all %in% labels)) {stop("'context.all' includes variable names not in 'labels'")}
      context.all <- which(labels %in% context.all)
    }
    
    if (is.numeric(context.all)) {
      if (!all(context.all %in% (1:p))) {stop("'context.all' contains elements that are smaller than 1 or larger than 'p'")}
      if (!all(tiers[context.all]==min(tiers))) {stop("'context.all' variables must be in the first tier")}
    } else {
      stop("'context.all' must be an integer vector or character vector")
    }
  }
  
  if (!is.null(context.tier)) {
    if (is.character(context.tier)) {
      if (!all(context.tier %in% labels)) {stop("'context.tier' includes variable names not in 'labels'")}
      context.tier <- which(labels %in% context.tier)
    }
    if (is.numeric(context.tier)) {
      if (!all(context.tier %in% 1:p)) {stop("'context.tier' contains elements that are smaller than 1 or larger than 'p'")}
    } else {
      stop("'context.tier' must be a numeric or character vector")
    }
  }
  
  if ( !is.null(context.tier) & !is.null(context.all) ) {
    if (length(intersect(context.tier, context.all)) > 0) {
      stop(paste("The following variables are in both 'context.tier' and 'context.all': ",
                 paste(intersect(context.tier, context.all), collapse=",")))
    }
  }
  
  
  ##Change1 ends
  
  
  
  
  
  
  ## Check that the type is a valid one
  type <- match.arg(type)
  if (type == "anytime" && m.max == Inf)
    stop("To use the Anytime FCI you must specify a finite 'm.max'.")
  if (type == "adaptive" && m.max != Inf)
    stop("To use the Adaptive Anytime FCI you must not specify 'm.max'.")
  
  if (conservative && maj.rule)
    stop("Choose either conservative FCI or majority rule FCI")
  
  
  
  
  #Change2 starts
  #Do we really need this step?
  if ((conservative) && (maj.rule))
    stop("Choose either conservative FCI or majority rule FCI")
  #Change2 ends
  
  
  
  
  ## Check that jci background knowledge is valid
  jci <- match.arg(jci)
  ## Check whether contextVars is valid
  if( !is.null(contextVars) && length(contextVars) > 0 ) {
    if( !is.numeric(contextVars) )
      stop("contextVars has to be a vector of integers in {1,2,..,p}, where p is the number of variables.")
    if( !all(sapply(contextVars, function(i) i == as.integer(i))) )
      stop("contextVars has to be a vector of integers in {1,2,..,p}, where p is the number of variables.")
    if( min(contextVars) < 1 || max(contextVars) > p )
      stop("contextVars has to be a vector of integers in {1,2,..,p}, where p is the number of variables.")
  }
  ## Set fixed edges from JCI assumption 3 if asked for
  if( jci == "123" && length(contextVars) > 0 ) {
    if( any(is.null(fixedEdges)) )
      fixedEdges<-matrix(FALSE,p,p)
    fixedEdges[contextVars,contextVars]<-TRUE
  }
  
  cl <- match.call()
  if (verbose) cat("Compute Skeleton\n================\n")
  
  skel <- tskeleton(suffStat, indepTest, alpha, labels, p,
                    method = c("stable", "original"), m.max = Inf,
                    fixedGaps = NULL, fixedEdges = NULL, NAdelete = TRUE,
                    tiers = NULL, verbose = FALSE)
  skel@call <- cl # so that makes it into result
  G <- as(skel@graph, "matrix")
  sepset <- skel@sepset
  pMax <- skel@pMax
  n.edgetestsSKEL <- skel@n.edgetests
  max.ordSKEL <- skel@max.ord
  allPdsep <- NA
  tripleList <- NULL
  
  if (doPdsep) {
    if (verbose) cat("\nCompute tPDSEP\n=============\n")
    #    pc.ci <- pc.cons.intern(skel, suffStat, indepTest,
    #                            alpha = alpha, version.unf = c(1,1),
    #                            maj.rule = FALSE, verbose = verbose)
    pc.ci <- list(unfTripl = c(), sk = skel)
    ## Recompute (sepsets, G, ...):
    pdsepRes <- tpdsep(skel@graph, suffStat, indepTest = indepTest, p = p,
                       sepset = pc.ci$sk@sepset, alpha = alpha, pMax = pMax,
                       m.max = if (type == "adaptive") max.ordSKEL else m.max,
                       pdsep.max = pdsep.max, NAdelete = NAdelete,
                       unfVect = pc.ci$unfTripl, # "tripleList.pdsep"
                       fixedEdges = fixedEdges, verbose = verbose, tiers = tiers)
    
    ## update the graph & sepset :
    G <- pdsepRes$G
    sepset <- pdsepRes$sepset
    pMax <- pdsepRes$pMax
    allPdsep <- pdsepRes$allPdsep
    n.edgetestsPD <- pdsepRes$n.edgetests
    max.ordPD <- pdsepRes$max.ord
    
    #--------code testing1 starts---------------#
    cat("show me the n.edgetestsPD from tpdsep...", n.edgetestsPD, "\n")
    #--------code testing1 starts---------------#
    
    if (conservative || maj.rule) {
      if (verbose)
        cat("\nCheck v-structures conservatively\n=================================\n")
      tmp.pdsep <- new("pcAlgo", graph = as(G, "graphNEL"), call = cl,
                       n = integer(0), max.ord = as.integer(max.ordSKEL),
                       n.edgetests = n.edgetestsSKEL, sepset = sepset,
                       pMax = pMax, zMin = matrix(NA, 1, 1))
      sk. <- tpc.cons.intern(tmp.pdsep, suffStat, indepTest, alpha,
                             verbose = verbose, version.unf = c(1, 1),
                             maj.rule = maj.rule)
      tripleList <- sk.$unfTripl
      ## update the sepsets
      sepset <- sk.$sk@sepset
    }
  }
  else {## !doPdsep : "do not Pdsep"
    n.edgetestsPD <- 0
    max.ordPD <- 0
    allPdsep <- vector("list", p)
    if (conservative || maj.rule) {
      if (verbose)
        cat("\nCheck v-structures conservatively\n=================================\n")
      nopdsep <- tpc.cons.intern(skel, suffStat, indepTest, alpha,
                                 verbose = verbose, version.unf = c(2, 1),
                                 maj.rule = maj.rule)
      tripleList <- nopdsep$unfTripl
      ## update the sepsets
      sepset <- nopdsep$sk@sepset
      
      #------------change4 starts
      ### INSERT 2 AT amat[j,i] WHEREVER THERE IS AN EDGE CONNECTING FUTURE at j TO PAST at i ###
      #storage.mode(G) <- "integer" # (TRUE, FALSE) -->  (1, 0)
      #for(i in 1:p) for(j in i:p) if(G[i,j]==1 && tiers[j] > tiers[i]) { 
      #  G[i,j] <- 2
      #  G[j,i] <- 1
      #}
      ##################
      #------------change4 ends
    }
  }
  if( !selectionBias )
    rules[5:7] <- FALSE
  else {
    if( jci != "0" )
      stop( 'The current JCI implementation does not support selection bias (use selectionBias=FALSE instead).' )
  }
  if (verbose)
    cat("\nDirect edges:\n=============\nUsing rules:", which(rules),
        "\nCompute collider:\n")
  
  #code test: check the G matrix for depicting:
  #cat("G IN The End of tfci algorithm", G, "\n")
  
  
  #------------change5 starts
  G1 <- G  #To keep G as a property, the G1 was created to be the proxy of G.
  ### INSERT 2 AT amat[j,i] WHEREVER THERE IS AN EDGE CONNECTING FUTURE at j TO PAST at i ###
  storage.mode(G1) <- "numeric" # (TRUE, FALSE) -->  (1, 0)
  for(i in 1:p) for(j in i:p) if(G1[i,j] == 1 && tiers[j] > tiers[i]) { 
    G1[i,j] <- 2
    G1[j,i] <- 1
  }
  ##################
  #------------change5 ends
  
  res <- udag2pag(pag = G1, sepset, rules = rules, unfVect = tripleList,
                  jci = jci, contextVars = contextVars, verbose = verbose)
  
  #code test: check the res matrix for depicting:
  #cat("res IN The End of tfci algorithm", res, "\n")
  
  
  
  colnames(res) <- rownames(res) <- labels
  
  
  new("fciAlgo", amat = res, call = cl, n = integer(0),
      max.ord = as.integer(max.ordSKEL),
      max.ordPDSEP = as.integer(max.ordPD),
      n.edgetests = n.edgetestsSKEL, n.edgetestsPDSEP = n.edgetestsPD,
      sepset = sepset, pMax = pMax, allPdsep = allPdsep)
  
  #code test: check the amat matrix for depicting:
  #cat("AMAT IN The End of tfci algorithm", res, "\n")
  
  
} ## {fci}