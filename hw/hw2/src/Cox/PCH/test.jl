include("PCH.jl")
include("../../est_survival.jl")

sym(M::Matrix{Float64}) = (M' + M) / 2

# Read data from R
using RCall

R"""
library(KMsurv)
library(survival)
library(rcommon)
data(larynx)
L <- model.matrix(time ~ as.factor(stage) -1 + age, data = larynx)[,-1]
colnames(L) <- c(paste0('stage',2:4),'age')
lt <- larynx$time
lv <- larynx$delta
""";
@rget L lt lv;

const t = lt
const v = lv * 1.
const X = L
const N = length(t)

### PCH
P = size(X,2)
#grid = [0; quantile(t,linspace(0,1,10))]
grid = collect(linspace(0,maximum(t),50))
J = length(grid) - 1
Σᵦ = sym(inv(X'X)) * 1.
priorβ = PCH.Priorβ(fill(0.,P), eye(P)*10., Σᵦ) 
priorλ = PCH.Priorλ(fill(.1,J), fill(.1,J), eye(Float64,J)*1E-1)
@time m2 = PCH.pch(t,X,v,grid,priorβ,priorλ,2000,10000,printFreq=500);
PCH.plot(m2,"beta",collect(1:P));
println(R"coxph(Surv(lt,lv)~L)")
s2 = PCH.summary(m2)
println(s2)

S = est_survival(m2, grid, 60., mean)
PCH.plotsurv(grid, S, lwd=3, col_l=["black","blue","red","green"],
            fg="grey",xlab="time", ylab="Survival Probability",
            main="Probability of Survival for different Stages",addlines=true)
