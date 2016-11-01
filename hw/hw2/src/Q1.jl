using RCall
include("Cox.jl")
Diag(v::Vector{Float64}) = convert(Matrix{Float64}, Diagonal(v))

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

B = 10000; burn = 20000; Σ = Diag([.005,.01,.01])#eye(3) * .005
@time m1 = Cox.coxph_weibull(t,x,d,Σ,B=B,burn=burn);
Cox.summary(m1)
Cox.plot(m1);

R"""
data(larynx)
L <- as.matrix(larynx)
"""
@rget L
@time m = Cox.coxph_weibull(L[:,2],L[:,[1,3,4]],L[:,5],
                            Diag([1E-2,1E-2,1E-2,1E-5,1E-5]),
                            B=10000,burn=30000);
Cox.summary(m)
println(R"coxph(Surv(time,delta) ~ ., data=larynx)")
Cox.plot(m);
