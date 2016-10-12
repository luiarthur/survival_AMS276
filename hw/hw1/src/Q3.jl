using RCall, MCMC, Distributions

R"""
library(KMsurv)   # To get the datasets in K-M
library(survival) # R functions
library(rcommon)  # devtools::install_github('luiarthur/rcommon')
data(tongue)
tongue$time <- tongue$time/10 # to obtain computational stability
""";

R"my.pairs(tongue)";
R"my.pairs(cbind(tongue[,1],log(tongue[,2]),tongue[,3]))";

aneuploid = R"tongue[which(tongue$type==1),]" 
diploid = R"tongue[which(tongue$type==2),]" 
