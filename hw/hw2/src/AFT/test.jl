include("../AFT/AFT.jl")
using RCall

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
""";
@rget tongue;

const t = tongue[:time];
const d = tongue[:delta];
const x = reshape(tongue[:type],length(t),1);
const N = length(t);

@time m2 = AFT.aft(t, [ones(N) x], d, [.1, .05], .1, B=2000, burn=10000, model="weibull");
b = hcat(map(m -> m.beta, m2)...)'; println(mean(b,1))
@rput b; R"plotPosts(b)";
println(R"summary(survreg(Surv(time,delta) ~ ., data=tongue))")

