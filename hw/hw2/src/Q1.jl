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
""";
@rget tongue;


### Analysis
const t = tongue[:time];
const d = tongue[:delta];
const x = reshape(tongue[:type],length(t),1);
const N = length(t);

srand(276);

### Cox
#=
Σ = Cox.Diag([.005,.01,.01])
@time m1 = Cox.coxph_weibull(t,x,d,Σ,B=2000,burn=20000);
println(Cox.Parametric.summary(m1))
Cox.Parametric.plot(m1);
println(R"coxph(Surv(time,delta) ~ type, data=tongue)")
=#

### PCH
J = 10
grid = [0; quantile(t,linspace(0,1,J))]
sym(X::Matrix{Float64}) = (X' + X) / 2
Σᵦ = sym(inv(x'x)) * 5
priorβ = Cox.Priorβ([0.], eye(1)*10., Σᵦ)
priorλ = Cox.Priorλ(fill(.1,J), fill(.1,J), eye(Float64,J)*1E-1)
@time m2 = Cox.pch(t,x,d,grid,priorβ,priorλ,2000,10000,printFreq=500)
s2 = Cox.PCH.summary(m2)

x0 = [[0.], [1.]]
mean_S_pch = Cox.PCH.est_survival(m2, grid, x0, mean)

Cox.PCH.plotsurv(grid, mean_S_pch, lwd=3, col_l=["blue","orange"],
             fg="grey",xlab="time", ylab="Survival Probability",
             addlines=true)
