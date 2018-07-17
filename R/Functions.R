#' SimCommunityAssembly
#'
#' @param sims: The number of simulations; sims
#' @param N: The size of the regional communities, can be a single number or a prior range
#' @param local: size of the local communities, can be a single number, a proportion of the regional pool (btw 0 and 1), or a prior range of proportions 
#' @param traitsim: either 'BM' or 'OU'
#' @param comsim: either "netural", "filtering", "competition"
#' @param lambda: speciation rate, can be single number or prior range, default is a unifrom prior between 0.5 and 2.0
#' @param eps: the extinction fraction, extinction (mu) = lambda * eps, default is a uniform prior between 0.2 and 0.8
#' @param sig2: rate of character change, can be single value or prior range, defualt is a unifom prior range between 1.0 and 10.0
#' @param alpha: pull to trait optimum, defualt is a uniform prior range between 0.01 and 0.2
#' @param tau: stregnth of the community assembly model, default is a uniform prior range between 5.0 and 5.0
#'
#' @return output is a list of two matrices, one containing all parameter values for each simulation, and the other containing all 
#'         summary statistics. In both matrices, each row corresponds to one simulation
#' @export
#'
#' @examples 
SimCommunityAssembly <- function(sims, N, local,
                                 traitsim,
                                 comsim,
                                 lambda = c(0.05, 2.0),
                                 eps = c(0.2, 0.8),
                                 sig2 = c(1, 10),
                                 alpha = c(0.01, 0.2),
                                 tau = c(5, 50)) {

  #check all parameters and reply with error message if any incomplete
  if (missing(sims))
    stop("'sims' is missing with no default", call. = F)
  if (missing(N))
    stop("'N' is missing with no default", call. = F)
  if (missing(local))
    stop("'local' is missing with no default", call. = F)
  if (missing(traitsim))
    stop("'traitsim' is missing with no default", call. = F)
  if (!any(traitsim == c("BM", "OU")))
    stop("Method must be 'BM' or 'OU'.",
         call. = F)
  if (missing(comsim))
    stop("'comsim' is missing with no default", call. = F)
  if (!any(comsim == c("neutral", "filtering", "competition")))
    stop("Method must be 'neutral', 'filtering', or 'competition'.",
         call. = F)

  ##empty vectors declared to hold info that will be output
  params <- matrix(NA, sims, 12)
  colnames(params) <- c("sim#",       #1
                        "comsim",     #2
                        "traitsim",   #3
                        "N",          #4
                        "local",      #5
                        "lambda",     #6
                        "mu",         #7
                        "sig2",       #8
                        "alpha",      #9
                        "tau",        #10
                        "avg.Pa",     #11
                        "avg.rej")    #12

  summary.stats <- matrix(NA, nrow = sims, 30)
  colnames(summary.stats) <-   c("Mean.BL", 	#1
                                 "Var.BL",		#2
                                 "Mean.Tr",		#3
                                 "Var.Tr",		#4
                                 "Moran.I",		#5
                                 "Age",			  #6
                                 "Colless",		#7
                                 "Sackin",		#8
                                 "nLTT",			#9
                                 "Msig",			#10
                                 "Cvar",			#11
                                 "Svar",			#12
                                 "Shgt",			#13
                                 "Dcdf",			#14
                                 "Skw",			  #15
                                 "Krts",			#16
                                 "MnRegTr",		#17
                                 "VrRegTr",		#18
                                 "MnTrDif",		#19
                                 "VrTrDif",		#20
                                 "MnRegBl",		#21
                                 "VarRegBl",	#22
                                 "MnBlDif",		#23
                                 "VarBlDif",	#24
                                 "Amp1",			#25
                                 "Amp2",			#26
                                 "BiCoef",		#27
                                 "BiRatio",		#28
                                 "Mode1",		  #29
                                 "Mode2")		  #30

  ##initiate progress bar
  pb <- txtProgressBar(min = 0, max = sims, style = 3)

  for (i in 1:sims){

    #drawn N if prior
    if (length(N) > 1)
      N <- runif(1, N[1], N[2])
    #draw lambda
    if (length(lambda) > 1)
      lambda.drawn <- runif(1, lambda[1], lambda[2])
    #draw eps
    if (length(eps) > 1)
      eps.drawn <- runif(1, eps[1], eps[2])
    #draw sig2
    if (length(sig2) > 1)
      sig2.drawn <- runif(1, sig2[1], sig2[2])
    #draw alpha
    if (length(alpha) > 1)
      alpha.drawn <- runif(1, alpha[1], alpha[2])
    #draw tau
    if (length(tau) > 1)
      tau.drawn <- runif(1, tau[1], tau[2])

    #Simulate Regional Tree
    if (is.integer(local)){
      n <- local
    }else{
      n <- N * local
    }
    mu <- lambda.drawn*eps.drawn
    regional.tree <- TreeSim::sim.bd.taxa(n=N, numbsim=1, lambda=lambda.drawn, mu=mu, complete=FALSE)[[1]]

    #Simulate Regional Traits
    if (traitsim == "BM") {
      traits <- geiger::sim.char(regional.tree, par=sig2.drawn, nsim=1, model="BM", root=0)[,,1]
    }else {
      traits <- phytools::fastBM(regional.tree, a=0, theta=0, alpha=alpha.drawn, sig2=sig2.drawn, internal=FALSE)
    }
    #store all regional traits here for calculating summary stats
    regional.traits <- traits


    #Simulate Community under given method
    rej = 0
    if (comsim == "neutral"){
      local.traits <- sample(traits, n)
      probs = c(1,1)
    }
    if (comsim == "filtering" || comsim == "competition"){
      #sample first species in communtiy; gets in no problem
      Xi <- sample(traits,1)
      local.traits <- Xi
      #remove this species/trait value from regional pool to be sample from
      traits <- traits[!traits==Xi]

      #Continuosly sample regional pool until local community is at required size as determined by 'local'
      probs = c()
      while (length(local.traits) < n && !rej > 5000) {
        Xj <- sample(traits,1)
        #calculate the probability of acceptance for this species into community, determined by Tau and the mean of the local trait values
        Pj <- (exp(-((mean(local.traits)-Xj)^2)/tau.drawn))

        if (comsim == "competition"){
          Pj <- 1 - Pj
        }

        if (Pj > runif(1,0,1)) {
          local.traits <- c(local.traits, Xj)
          traits <- traits[!traits==Xj]
          probs = c(probs, Pj)
        }else{
          rej=rej+1
          probs = c(probs, Pj) }
      }
    }

    #construct local tree
    local.tree <- ape::drop.tip(regional.tree,setdiff(regional.tree$tip.label, names(local.traits)))
    #local.tree <- keep.tip(regional.tree, names(local.traits))

    #Store all parameters for each simulation and data produced
    params[i, ] <- c(i, comsim, traitsim, N, n, lambda.drawn, mu, sig2.drawn, alpha.drawn, tau.drawn, mean(probs), rej)

    ##calculate all summary statistics
    summary.stats[i, ] <- CalcSummaryStats(regional.tree, local.tree, regional.traits, local.traits)

    Sys.sleep(.1)
    # update progress bar
    setTxtProgressBar(pb, i)
  }

  output <- list(data.frame(params), summary.stats)
  names(output) <- c("params", "summary.stats")
  close(pb)
  return(output)

}

#' CalcSummaryStats
#'
#' @param regional.tree: regional community phylogenetic tree, of class "phylo"
#' @param local.tree: local community phylogenetic tree, of class "phylo"
#' @param regional.traits: trait values for each species in regional community, vector or list
#' @param local.traits: trait values for each species in local community, vector or list
#'
#' @return vector of summary statistics
#' @export
#'
#' @examples: 
CalcSummaryStats <- function(regional.tree,
                             local.tree,
                             regional.traits,
                             local.traits){

  if (class(regional.tree)!="phylo")
    stop("'regional.tree' is not a phylo object", call. = F)
  if (class(local.tree)!="phylo")
    stop("'regional.tree' is not a phylo object", call. = F)
  if (length(regional.traits)!=length(regional.tree$tip.label))
    stop("Number of regional traits does not match the number of tips on the regional tree", call. = F)
  if (length(local.traits)!=length(local.tree$tip.label))
    stop("Number of local traits does not match the number of tips on the local tree", call. = F)

  ## 1. Mean branch length in local tree
  mean.branch.local.tree <- mean(local.tree$edge.length)

  ## 2. Variance of branch lengths in local tree
  var.branch.local.tree <- var(local.tree$edge.length)

  ## 3. Mean trait value in local community
  mean.trait.local <- mean(local.traits)

  ## 4. Variance of trait values in local community
  var.trait.local <- var(local.traits)

  ## 5. Moran's autocorrelation index
  w <- 1/cophenetic(local.tree)
  diag(w) <- 0
  MI <-  ape::Moran.I(local.traits, w)$observed

  ## 6. Maximum node depth
  max.depth <- max(ape::node.depth.edgelength(local.tree))

  ## 7. Colless index
  colless <- phyloTop::colless.phylo(local.tree)

  ## 8. Sackin index
  sackin <- phyloTop::sackin.phylo(local.tree)

  ## 9. lineage through time plot differences between regional and local tree
  nLTT <- nLTT::nLTTstat_exact(regional.tree, local.tree)

  ## 10. Mean of squared phylogenetic independent contrasts in local communtiy
  PIC <- ape::pic(local.traits, local.tree, var.contrasts=T)
  Msig <- mean(PIC[,1]^2)

  ## 11. Standard devation of PICs in local over the mean PICs in local
  Cvar <- sd(PIC[,1])/mean(PIC[,1])

  ## 12. slope of linear model between abs(PICs) and the expected variance in PICs
  Svar <- lm(abs(PIC[,1]) ~ PIC[,2])$coefficients[2]

  ## 13. Slope of linear model between the absolute magnitude of the standardized independent contrasts
  ##     and the height above the root of the node at which they were being compared
  Shgt.lm <- geiger::nh.test(local.tree, local.traits[local.tree$tip.label], regression.type="lm", show.plot=FALSE)
  Shgt <- Shgt.lm$coefficients[2]

  ## 14. D statistic from KS.test between PIC of local communtiy and expect variance under BM
  expCont.BM <- rnorm(n=length(local.traits), mean=0, sd=sqrt(mean(PIC[,1]^2)))
  Dcdf <- ks.test(PIC[,1], expCont.BM)$statistic

  ## 15. skeness of local community traits
  skew <- modes::skewness(local.traits)

  ## 16. Kurtosis of local community traits
  kurt <- modes::kurtosis(local.traits, finite=T)

  ## 17. Mean of regional trait values
  mean.trait.reg <- mean(regional.traits)

  ## 18. Variance of regional trait values
  var.trait.reg <- var(regional.traits)

  ## 19. Difference in mean of trait values between regional and local communities
  mean.trait.dif <- mean.trait.reg - mean.trait.local

  ## 20. Difference in variance of trait values between regional and local communities
  var.trait.dif <- var.trait.reg - var.trait.local

  ## 21. Mean branch length in regional tree
  mean.branch.reg.tree <- mean(regional.tree$edge.length)

  ## 22. Variance of branch lengths in regional tree
  var.branch.reg.tree <- var(regional.tree$edge.length)

  ## 23. Difference between mean of branch lengths between local and regional tree
  mean.branch.dif <- mean.branch.reg.tree - mean.branch.local.tree

  ## 24. Difference between variance of branch lengths between local and regional tree
  var.branch.dif <- var.branch.reg.tree - var.branch.local.tree

  ## 25. density peak location in local community traits
  amp <- modes::amps(local.traits)$Peaks
  amp1 <- amp[1]

  ## 26. density peak amplitude in local community traits
  amp2 <- amp[2]

  ## 27. Bimodality coefficient
  bimo.coef <- modes::bimodality_coefficient(local.traits)

  ## 28. Bimodality ratio
  bimodal <- modes::bimodality_ratio(local.traits)
  if (is.na(bimodal)){
    bimodal <- 0
  }

  ## 29.
  modes <- modes::modes(local.traits)
  modes1 <- modes[1]

  ## 30.
  modes2 <- modes[2]

  stats <-  c(mean.branch.local.tree, 	#1
              var.branch.local.tree,		#2
              mean.trait.local,		#3
              var.trait.local,		#4
              MI,		#5
              max.depth,			  #6
              colless,		#7
              sackin,		#8
              nLTT,			#9
              Msig,			#10
              Cvar,			#11
              Svar,			#12
              Shgt,			#13
              Dcdf,			#14
              skew,			  #15
              kurt,			#16
              mean.trait.reg,		#17
              var.trait.reg,		#18
              mean.trait.dif,		#19
              var.trait.dif,		#20
              mean.branch.reg.tree,		#21
              var.branch.reg.tree,	#22
              mean.branch.dif,		#23
              var.branch.dif,	#24
              amp1,			#25
              amp2,			#26
              bimo.coef,		#27
              bimodal,		#28
              modes1,		  #29
              modes2)		  #30

  names(stats) <- c("Mean.BL", 	#1
                    "Var.BL",		#2
                    "Mean.Tr",		#3
                    "Var.Tr",		#4
                    "Moran.I",		#5
                    "Age",			  #6
                    "Colless",		#7
                    "Sackin",		#8
                    "nLTT",			#9
                    "Msig",			#10
                    "Cvar",			#11
                    "Svar",			#12
                    "Shgt",			#13
                    "Dcdf",			#14
                    "Skw",			  #15
                    "Krts",			#16
                    "MnRegTr",		#17
                    "VrRegTr",		#18
                    "MnTrDif",		#19
                    "VrTrDif",		#20
                    "MnRegBl",		#21
                    "VarRegBl",	#22
                    "MnBlDif",		#23
                    "VarBlDif",	#24
                    "Amp1",			#25
                    "Amp2",			#26
                    "BiCoef",		#27
                    "BiRatio",		#28
                    "Mode1",		  #29
                    "Mode2")		  #30

  return(stats)

}


