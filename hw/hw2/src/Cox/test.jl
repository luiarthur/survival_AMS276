include("Cox.jl")
include("../../../hw1/src/AFT.jl")
using RCall

R"""
library(KMsurv)
library(survival)
data(larynx)
L <- as.matrix(larynx)
"""
@rget L;

const t = L[:,2]
const X = L[:,[1,3,4]]
const v = L[:,5]

@time m = Cox.coxph_weibull(t,X,v,
                            Cox.Diag([1E-2,1E-2,1E-2,1E-5,1E-5]),
                            B=10000,burn=30000);
println(R"summary(survreg(Surv(time,delta) ~ ., data=larynx))")
println(R"summary(coxph(Surv(time,delta) ~ ., data=larynx))")
s = Cox.summary(m)
Cox.plot(m);

#=
update AFT.aft so that it can have matrix proposal
=#


@time m2 = AFT.aft(t, X, v, css=.5, csb=.05, B=2000, burn=100000);
b = hcat(map(m -> m.beta, m2)...)';
mean(b,1)
@rput b;
R"plotPosts(b)";

