include("Cox.jl")
using RCall

R"""
library(KMsurv)
library(survival)
data(larynx)
L <- as.matrix(larynx)
"""
@rget L;
@time m = Cox.coxph_weibull(L[:,2],L[:,[1,3,4]],L[:,5],
                            Cox.Diag([1E-2,1E-2,1E-2,1E-5,1E-5]),
                            B=10000,burn=30000);
Cox.summary(m)
println(R"coxph(Surv(time,delta) ~ ., data=larynx)")
Cox.plot(m);
