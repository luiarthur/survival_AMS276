include("Cox.jl")
include("../AFT/AFT.jl")
using RCall

R"""
library(KMsurv)
library(survival)
library(rcommon)
data(larynx)
L <- model.matrix(time ~ as.factor(stage) -1 + age, data = larynx)[,-1]
lt <- larynx$time
lv <- larynx$delta
""";
@rget L lt lv;

const t = lt
const v = lv * 1.
const X = L
const N = length(t)

@time m = Cox.coxph_weibull(t,X,v,
                            Cox.Diag([.1,.1,.1,1,.1,1]*1E-4),
                            B=20000,burn=100000);
s = Cox.summary(m)
Cox.plot(m);
println(R"summary(coxph(Surv(larynx$time,larynx$delta) ~ L,method='breslow'))")


@time m2 = AFT.aft(t, [ones(N) X], v, [.0001,.0001,.0001,.0001,.0001], .3, B=20000, burn=300000, model="weibull");
b = hcat(map(m -> m.beta, m2)...)'; mean(b,1)
s = hcat(map(m -> m.sig, m2)...)'; mean(s)
@rput b; R"plotPosts(b)";
println(R"summary(survreg(Surv(larynx$time,larynx$delta) ~ L))")
