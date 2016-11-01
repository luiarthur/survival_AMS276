using RCall
include("Cox.jl")

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
print(coxph(Surv(time,delta) ~ type, data=tongue))
"""
@rget tongue
println()

const t = tongue[:time]
const d = tongue[:delta]
const x = reshape(tongue[:type],length(t),1)

srand(276);

B = 10000; burn = 20000; Σ = eye(3) * .005
@time m1 = Cox.coxph_weibull(t,x,d,Σ,B=B,burn=burn);
s1 = Cox.summary_cox(m1)
println(s1)

beta = map(p -> p.β[1], m1.params)
alpha = map(p -> p.α, m1.params)
lambda = map(p -> p.λ, m1.params)

@rput beta lambda alpha

R"plotPosts(cbind(beta,alpha,lambda),show.x=F)";
