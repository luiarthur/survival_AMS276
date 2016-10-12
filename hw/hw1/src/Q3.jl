using RCall, MCMC, Distributions

R"""
library(KMsurv)   # To get the datasets in K-M
library(survival) # R functions
library(rcommon)  # devtools::install_github('luiarthur/rcommon')
data(tongue)
tongue$time <- tongue$time/10 # to obtain computational stability
""";

R"my.pairs(tongue)";

aneuploid = R"tongue[which(tongue$type==1),]" 
diploid = R"tongue[which(tongue$type==2),]" 
