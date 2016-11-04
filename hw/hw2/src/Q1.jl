include("Cox/Cox.jl")

# Read data from R
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

"""
@rget tongue


### Analysis
const t = tongue[:time]
const d = tongue[:delta]
const x = reshape(tongue[:type],length(t),1)

srand(276);

B = 10000; burn = 20000; Σ = Cox.Diag([.005,.01,.01])
@time m1 = Cox.coxph_weibull(t,x,d,Σ,B=B,burn=burn);
Cox.summary(m1)
Cox.plot(m1);

println(R"coxph(Surv(time,delta) ~ type, data=tongue)")
println(R"summary(survreg(Surv(time,delta) ~ type, dist='weibull', data=tongue))")
