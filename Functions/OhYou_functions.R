# Functions for running OU paper simulations

# Simulating trees
# Using diversitree function: trees

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

# Yule tree simulations
yule.trees <- function(ntaxa) {
  simulate.trees(lambda = 1, mu = 0, ntaxa)
}

# Birth death tree simulations with low mu
bdlow.trees <- function(ntaxa) {
  simulate.trees(lambda = 1, mu = 0.25, ntaxa)
}

# Birth death tree simulations with medium mu
bdmid.trees <- function(ntaxa) {
  simulate.trees(lambda = 1, mu = 0.5, ntaxa)
}

# Birth death tree simulations with high mu
bdhigh.trees <- function(ntaxa) {
  simulate.trees(lambda = 1, mu = 0.75, ntaxa)
}

# Running all tree simulations
all.trees <- function(ntaxa) {
  yule.trees(ntaxa)
  bdlow.trees(ntaxa)
  bdmid.trees(ntaxa)
  bdhigh.trees(ntaxa)
}
