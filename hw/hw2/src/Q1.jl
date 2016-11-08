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

### Cox weibull
Σ = Cox.Diag([.005,.01,.01])
@time m1 = Cox.coxph_weibull(t,x,d,Σ,B=2000,burn=20000);
s1 = Cox.Parametric.summary(m1)
println(s1)
Cox.Parametric.plot(m1);
t0 = collect(linspace(0,maximum(t),1000))
mean_S_weib = Cox.Parametric.est_survival(m1, t0, [[0.],[1.]], mean)
Cox.Parametric.plotsurv(t0,mean_S_weib, lwd=3, col_l=["blue","orange"],
                        fg="grey",xlab="time", ylab="Survival Probability")
R"lines(survfit(Surv(time,delta) ~ type, data = tongue))"

### PCH
J = 10
grid = [0; quantile(t,linspace(0,1,J))]
sym(X::Matrix{Float64}) = (X' + X) / 2
Σᵦ = sym(inv(x'x)) * 30
priorβ = Cox.Priorβ([0.], eye(1)*10., Σᵦ)
priorλ = Cox.Priorλ(fill(.1,J), fill(.1,J), eye(Float64,J)*1E-1)
@time m2 = Cox.pch(t,x,d,grid,priorβ,priorλ,2000,10000,printFreq=500);
s2 = Cox.PCH.summary(m2)
println(s2)
Cox.PCH.plot(m2,"beta", [1]);

x0 = [[0.], [1.]]
mean_S_pch = Cox.PCH.est_survival(m2, grid, x0, mean)

Cox.PCH.plotsurv(grid, mean_S_pch, lwd=3, col_l=["blue","orange"],
             fg="grey",xlab="time", ylab="Survival Probability",
             addlines=true)
R"lines(survfit(Surv(time,delta) ~ type, data = tongue))";

### GP
@time m3 = Cox.GammaProcess.gp(t,x,d,5.,.07,2000,1000, printFreq=500);
s3 = Cox.GammaProcess.summary(m3)
println(s3)
Cox.GammaProcess.plot(m3,"beta", [1]);
Cox.GammaProcess.plot(m3,"h",[1,2,3,4,5]);
Cox.GammaProcess.plot(m3,"h",[6,7,8,9,10]);

grid_gp = sort(unique([0;t]))
x0 = [[0.], [1.]]
mean_S_gp = Cox.GammaProcess.est_survival(m3, grid_gp, x0, mean)
Cox.PCH.plotsurv(grid_gp, mean_S_gp, lwd=3, col_l=["blue","orange"],
                 fg="grey",xlab="time", ylab="Survival Probability",
                 main="Probability of Survival for different Stages",addlines=true);
R"lines(survfit(Surv(time,delta) ~ type, data = tongue))";

###
println(R"coxph(Surv(time,delta) ~ type, data=tongue)")

R"par(mfrow=c(1,3))";
Cox.Parametric.plotsurv(t0,mean_S_weib, lwd=3, col_l=["blue","orange"],
                        fg="grey",xlab="time", ylab="Survival Probability");
R"title(main='Parametric Cox Model',cex.main=2)"
R"lines(survfit(Surv(time,delta) ~ type, data = tongue))";
R"""
legend("topright",legend=c("aneuplod","diploid"),col=c("orange","blue"),bty="n",
       lwd=5,cex=3,text.col='grey')
"""

Cox.PCH.plotsurv(grid, mean_S_pch, lwd=3, col_l=["blue","orange"],
                 fg="grey",xlab="time", ylab="Survival Probability",
                 addlines=true);
R"title(main='PCH Cox Model',cex.main=2)"
R"lines(survfit(Surv(time,delta) ~ type, data = tongue))";

Cox.PCH.plotsurv(grid_gp, mean_S_gp, lwd=3, col_l=["blue","orange"],
             fg="grey",xlab="time", ylab="Survival Probability",
             addlines=true);
R"title(main='GP Cox Model',cex.main=2)"
R"lines(survfit(Surv(time,delta) ~ type, data = tongue))";
R"par(mfrow=c(1,1))";

