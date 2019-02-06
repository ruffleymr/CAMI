# CAMI
### Community Assembly Model Inference

  For full Installation instructions and initial uses, refer to vignette ‘CAMI_Installation.Rmd’. Briefly, you can install the R package using devtools::install_github("ruffleymr/CAMI")

#### **Overview**

  CAMI employs a stochastic algorithm to simulate communities assembled under environmental filtering, competitive exclusion, and neutral species assembly processes. The algorithm was adapted from coevolutionary phenotypic matching and repulsion models. CAMI simultaneously considers the phylogenetic and phenotypic information from species in the local and regional communities and parameterizes the relative strength of the assembly processes to mimic strong to mild non-neutral assembly. CAMI implements a model-based inference procedure by using two approximate approaches, random forests (RF; Breiman 2001; Breiman & Cutler 2007) and approximate Bayesian computation (ABC; Csilléry et al. 2010). Additionally, because the strength of non-neutral assembly models is parameterized, the strength can be estimated using ABC.

#### **Data Simulation**

  For a single simulation of community assembly, first, a regional community phylogeny is simulated under a constant birth-death process with speciation and extinction parameters, until the desired number of regional species, N, is reached (Stadler 2011). Traits are evolved on the regional phylogeny (Revell 2012) under either a Brownian Motion (BM) or Ornstein-Uhlenbeck (OU) model of trait evolution (Butler & King 2004) characterized by the rate of character change and for OU models, the “strength of pull” to the trait optimum. Once the regional community exits with phylogenetic relationships and trait information, the assembly of the local community can begin.

  The assembly process uses the probabilities of species persisting in local communities for environmental filtering and competitive exclusion, and a rejection algorithm to stochastically assemble the local community. When a species colonizes the community, the probability of persistence is calculated, and the species is included in the local community if that probability is greater than a uniform random number between 0 and 1. Otherwise, the species is rejected from being in the local community. When a species is rejected from entering the community, it remains in the regional pool and is still able to colonize the local community again. In this case, the probability of persistence is recalculated, and the species has another chance to pass the rejection algorithm. As in the neutral model, the assembly process ends when the local community has reached capacity.

  All parameters mentioned are either fixed or drawn from a prior distribution. Information regarding the default prior distributions or fixed values for each parameter can be found in the help documentation for the R package ‘CAMI’ (github.com/ruffleymr/CAMI).

#### **Inference Procedure**

  For a single simulation of community assembly, a regional and local phylogeny and a regional and local distribution of trait values is returned. This information is summarized in 30 different summary statistics that capture information about the phylogeny, trait distributions, and phylogenetic signal within the traits of the local community (Komsta & Novomestky 2015, Janzen et al. 2015; Pennell et al. 2015; Deevi et al. 2016, Kendall et al. 2018, Paradis & Schliep 2018; see function CalcSummaryStats()). These summary statistics are then used for model selection and parameter estimation.

  To predict model probabilities from empirical data, we used two model selection approaches. The first approach uses a machine learning classification algorithm, random forests (RF; Breiman 1999; Liaw & Wiener 2002) to build a ‘forest’ of classification trees using the simulated summary statistics as predictor variables and the community assembly models as response variables. As a classifier is being built, RF is simultaneously measuring the ‘Out of Bag’ (OoB) error rates of the classifier by cross-validating each classification tree with a subset of the original data that was not used to make the tree in question. The OoB error rates measure how often the data are incorrectly classified. Additionally, RF quantifies the effect of including each summary statistic on the accuracy of the classifier through two variable importance measures, Mean Decrease in Accuracy (MDA) and Mean decrease in Gini Index (GINI) (Breiman 2002).

  RF is generally robust to noisy and/or overpowering predictor variables because each tree in the forest is constructed with only a subset of the data and multiple predictor variables are used at each node (Breiman & Cutler 2007). Our second approach, ABC, relies on the Euclidean distance between observed and simulated summary statistics to accept simulations into the posterior probability distribution of the models given the data (Csilléry et al. 2010). The support for each model then comes from the proportion of simulations from each model accepted into the posterior probability distribution. If there are summary statistics included that add a lot of noise to the classification process, ABC will lose power in distinguishing support between models. RF is able to measure which summary statistics are the most influential in distinguishing between the models, through importance measures such as MDA and GINI, thus we used this information to select a subset of 10 summary statistics to use in ABC model selection. ABC then predicts model probabilities using those statistics, a rejection algorithm, and a tolerance of 0.001 (Csilléry et al. 2012). The performance of ABC in classifying the data can be measured using a leave-one-out cross validation approach for model selection which results in model misclassification rates for each model.

#### **References**

Breiman, L. (2001). Random forests. Mach. Learn.45, 5-32

Breiman, L. (2002). Manual on setting up, using, and understanding random forests v3. 1. Stat. Dep. Univ. Calif. Berkeley, CA, USA.

Breiman, L. & Cutler, A. (2007). Random forests — Classification description: Random forests. http//stat-www.berkeley.edu/users/breiman/RandomForests/cc_home. htm.

Butler, M.A. & King, A.A. (2004). Phylogenetic Comparative Analysis: A Modeling Approach for Adaptive Evolution. Am. Nat.,164, 683-695.

Csilléry, K., Blum, M.G.B., Gaggiotti, O.E. & François, O. (2010). Approximate Bayesian Computation (ABC) in practice. Trends Ecol. Evol., 25, 410-418.

Csilléry, K., François, O. & Blum, M.G.B. (2012). Abc: An R package for approximate Bayesian computation (ABC). Methods Ecol. Evol., 3, 475–479.

Janzen, T., Höhna, S. & Etienne, R.S. (2015). Approximate Bayesian Computation of diversification rates from molecular phylogenies: Introducing a new efficient summary statistic, the nLTT. Methods Ecol. Evol., 6, 566–575.

Liaw, A. & Wiener, M. (2002). Classification and Regression by randomForest. R news., 2/3, 18-22.
Paradis, E. & Schliep, K. (2018). ape 5.0: an environment for modern phylogenetics and evolutionary analyses in R. Bioinformatics., 35, 526-528.

Pennell, M.W., FitzJohn, R.G., Cornwell, W.K. & Harmon, L.J. (2015). Model Adequacy and the Macroevolution of Angiosperm Functional Traits. Am. Nat., 186, E33–E50.
Revell, L.J. (2012). phytools: An R package for phylogenetic comparative biology (and other things). Methods Ecol. Evol., 3, 217–223.

Stadler, T. (2011). Simulating trees with a fixed number of extant species. Syst. Biol., 60, 676–684.
