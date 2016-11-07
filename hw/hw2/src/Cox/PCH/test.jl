include("PCH.jl")

# Read data from R
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

### PCH
J = 3
const P = size(X,2)
grid = [0; quantile(t,linspace(0,1,J))]
priorβ = PCH.Priorβ(fill(0.,P), eye(P)*10., eye(P)*1E-5)
priorλ = PCH.Priorλ(fill(.1,J), fill(.1,J), eye(Float64,J)*1E-5)
@time m2 = PCH.pch(t,X,v,grid,priorβ,priorλ,10000,1000,printFreq=100)
s2 = PCH.summary(m2)
PCH.plot(m2)

println(R"coxph(Surv(lt,lv) ~ L)")
