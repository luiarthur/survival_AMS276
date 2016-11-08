include("PCH.jl")
#include("../../est_survival.jl")

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
#grid = collect(linspace(0,maximum(t),10))
grid = [0; quantile(t, linspace(0,1,10))]
J = length(grid) - 1
Σᵦ = sym(inv(X'X)) * 1.5
priorβ = PCH.Priorβ(fill(0.,P), eye(P)*10., Σᵦ) 
priorλ = PCH.Priorλ(fill(.1,J), fill(.1,J), eye(Float64,J)*1E-1)
@time m2 = PCH.pch(t,X,v,grid,priorβ,priorλ,2000,30000,printFreq=500);
PCH.plot(m2,"beta",collect(1:P));
PCH.plot(m2,"lambda",collect(1:4));
println(R"coxph(Surv(lt,lv)~L)")
s2 = PCH.summary(m2)
println(s2)

x0 = [[0,0,0,60.],[1,0,0,60.],[0,1,0,60.],[0,0,1,60.]]
mean_S = PCH.est_survival(m2, grid, x0, mean)
q_025_S = PCH.est_survival(m2, grid, x0, s->quantile(s,.025))
q_975_S = PCH.est_survival(m2, grid, x0, s->quantile(s,.975))
PCH.plotsurv(grid, mean_S, lwd=3, col_l=["black","blue","red","green"],
             fg="grey",xlab="time", ylab="Survival Probability",
             main="Probability of Survival for different Stages",addlines=true)
PCH.plotsurv(grid, q_025_S, lwd=3, col_l=["black","blue","red","green"],
             fg="grey",xlab="time", ylab="Survival Probability",add=false,
             main="Probability of Survival for different Stages",addlines=true)
PCH.plotsurv(grid, q_975_S, lwd=3, col_l=["black","blue","red","green"],
             fg="grey",xlab="time", ylab="Survival Probability",add=false,
             main="Probability of Survival for different Stages",addlines=true)
