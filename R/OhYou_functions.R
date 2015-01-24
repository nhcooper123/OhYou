# Functions for running OU paper simulations

#----------------------------------
# Simulating trees
# Using diversitree function: trees
#----------------------------------

# Rename tip labels so they match other functions
rename.tips <- function(phy) {
  phy$tip.label <- paste("t", 1:Ntip(phy), sep = "")
  return(phy)
}

# Rejection algorithm to avoid simulating trees with
# fewer taxa than required
check.tree.size <- function(phy, lambda, mu, ntaxa) {
  while(length(phy$tip.label) != ntaxa) {
  	simulate.trees(lambda, mu, ntaxa)
  }  
}

# General tree simulation function
simulate.trees <- function(lambda, mu, ntaxa) {
  tree.sim <- trees(pars = c(lambda, mu), type = "bd", n = 1,
                    max.taxa = ntaxa, max.t = 10)
  tree.sim <- rename.tips(tree.sim[[1]])
  check.tree.size(tree.sim, lambda, mu, ntaxa)
  return(tree.sim)
}

# Write trees to file for later use
write.trees <- function(phy, ntaxa, treetype) {
  write.tree(phy, file = paste(treetype, ntaxa, ".tre", sep = ""), 
             append = TRUE)
}

# Yule tree simulations
yule.trees <- function(ntaxa, write.tree = FALSE) {
  tree.sim <- simulate.trees(lambda = 1, mu = 0, ntaxa)
  if(write.tree = TRUE) {
    write.trees(tree.sim, ntaxa, "yule")
  }
}

# Birth death tree simulations with low mu
bdlow.trees <- function(ntaxa, write.tree = FALSE) {
  simulate.trees(lambda = 1, mu = 0.25, ntaxa)
  if(write.tree = TRUE) {
    write.trees(tree.sim, ntaxa, "bdlow")
  }
}

# Birth death tree simulations with medium mu
bdmid.trees <- function(ntaxa, write.tree = FALSE) {
  simulate.trees(lambda = 1, mu = 0.5, ntaxa)
  if(write.tree = TRUE) {
    write.trees(tree.sim, ntaxa, "bdmid")
  }
}

# Birth death tree simulations with high mu
bdhigh.trees <- function(ntaxa, write.tree = FALSE) {
  simulate.trees(lambda = 1, mu = 0.75, ntaxa)
  if(write.tree = TRUE) {
    write.trees(tree.sim, ntaxa, "bdhigh")
  }
}

#--------------------------------------------------
# Adding error to branches leading to tips of trees
#--------------------------------------------------

id.tip.branches <- function(phy) {
  match(1:Ntip(phy), phy$edge[,2])
}

get.error <- function(phy, error) {
  error * max(node.depth.edgelength(phy))
}

get.tree.with.error <- function(phy, error) {
  phy.error <- phy
  phy.error$edge.length[id.tip.branches(phy.error)] <- phy.error$edge.length[id.tip.branches(phy.error)] 
                                                        + get.error(phy, error)
  return(phy.error)
}

#-----------------
# Simulating data
#-----------------

# Simulating data under Brownian model
bm.data <- function(phy, sigma) {
  rTraitCont(phy, model = "BM", sigma = sigma)
}

# Simulating data under OU model
# For now theta is hard coded as 0, same as the root state
# This matches the GEIGER model which does not estimate theta
# and instead uses the same theta as the root
ou.data <- function(phy, sigma, alpha) {
  rTraitCont(phy, model = "OU", sigma = sigma, alpha = alpha, theta = 0)
}

#------------------------------------------------------
# Fitting BM and OU models to simulated trees and data
# ML and Bayesian models if possible
#------------------------------------------------------


