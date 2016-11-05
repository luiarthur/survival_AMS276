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


@time m2 = AFT.aft(t, X, v, [.3,.15,.005,.005], .5, B=10000, burn=500000);
b = hcat(map(m -> m.beta, m2)...)'; mean(b,1)
@rput b; R"plotPosts(b)";
println(R"summary(survreg(Surv(time,delta) ~ ., data=larynx))")
