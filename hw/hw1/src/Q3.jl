using RCall, MCMC

R"""
library(KMsurv)   # To get the datasets in K-M"
library(survival) # R functions"
library(rcommon)  # devtools::install_github('luiarthur/rcommon')
tongue$time <- tongue$time/10 # to obtain computational stability
data(tongue)
""";

R"my.pairs(tongue)";

