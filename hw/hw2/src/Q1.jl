using RCall
include("../../hw1/src/AFT.jl")

R"""
library(KMsurv)
library(survival)
library(rcommon)
data(tongue)

tongue$time <- tongue$time/10 # for computational stability

tongue$delta <- ifelse(tongue$delta==1,1,0)
# 1: observed
# 2: right-censored

tongue$type <- ifelse(tongue$type==1,1,0)
# aneuplod = 1
#  diploid = 0

# see the result
survreg(Surv(time,delta) ~ type, dist='weibull', data=tongue)
coxph(Surv(time,delta) ~ type, data=tongue)
"""

@rget tongue

const t = tongue[:time]
const d = tongue[:delta]
const x = reshape(tongue[:type],length(t),1)
const B = 10000
const burn = 5000

srand(276);

@time m1 = AFT.aft(t, x, d, B=B, burn=burn);

b_post = hcat(map(x -> x.beta, m1)...)'
s_post = map(x -> x.sig, m1)

@rput b_post s_post;
R"plotPosts(b_post)";
R"plotPost(s_post)";
