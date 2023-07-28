# temporal-fast-casual-inference-research-project

Authors: Qixiang Chen and Daniel Malinsky, Columbia University.

This is an R implementation of the FCI algorithm with tiered background knowledge.

It is based heavily on the existing code from R packages 'pcalg' and 'tpc' (https://github.com/bips-hb/tpc).

See also implementation of FCI with tiered background knowledge in TETRAD (written in Java: https://github.com/cmu-phil/tetrad).

Syntax: 
> source("tfci.R")
> tfci(suffStat, indepTest, alpha = .05, p = 4, verbose=FALSE, tiers = c(1, 1, 1, 2)) ## see test examples

** Note: tiers must be specified in increasing order (so input variables in data matrix must be specified in increasing tier order).

Still in development!
