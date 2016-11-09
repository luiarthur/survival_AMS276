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
rgb=R"rgb"


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
t0 = collect(linspace(0,maximum(t),100))
x0 = [[0.], [1.]]
mean_S_weib = Cox.Parametric.est_survival(m1, t0, x0, mean)
param_025 = Cox.Parametric.est_survival(m1, t0, x0, m->quantile(m,.025))
param_975 = Cox.Parametric.est_survival(m1, t0, x0, m->quantile(m,.975))
Cox.PCH.plotCI(t0,[param_025[:,2] param_975[:,2]],col_area="orange",xlab="months")
Cox.PCH.plotCI(t0,[param_025[:,1] param_975[:,1]],col_area=rgb(0,0,1,.5),add=true)
Cox.PCH.plotsurv(t0, mean_S_weib, lwd=3, col_l=["blue","yellow"],
                 fg="grey",xlab="months", ylab="Survival Probability",
                 add=true,addlines=true);
R"lines(survfit(Surv(time,delta) ~ type, data = tongue),col='grey30',lwd=2)";

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

mean_S_pch = Cox.PCH.est_survival(m2, grid, x0, mean)
pch_025 = Cox.PCH.est_survival(m2, grid, x0, m->quantile(m,.025))
pch_975 = Cox.PCH.est_survival(m2, grid, x0, m->quantile(m,.975))
Cox.PCH.plotCI(grid, [pch_025[:,2] pch_975[:,2]],col_area="orange",xlab="months")
Cox.PCH.plotCI(grid, [pch_025[:,1] pch_975[:,1]],col_area=rgb(0,0,1,.5),add=true)
Cox.PCH.plotsurv(grid, mean_S_pch, lwd=3, col_l=["blue","yellow"],
                 fg="grey",xlab="months", ylab="Survival Probability",
                 add=true,addlines=true);
R"lines(survfit(Surv(time,delta) ~ type, data = tongue),col='grey',lwd=2)";


### GP
@time m3=Cox.GammaProcess.gp(t,x,d,10.,.07,2000,1000,
                             c=1.,η=1.,κ=1.,printFreq=500);
s3 = Cox.GammaProcess.summary(m3)
println(s3)
Cox.GammaProcess.plot(m3,"beta", [1]);
Cox.GammaProcess.plot(m3,"h",[1,2,3,4,5]);
Cox.GammaProcess.plot(m3,"h",[6,7,8,9,10]);

grid_gp = sort(unique([0;t]))
mean_S_gp = Cox.GammaProcess.est_survival(m3, grid_gp, x0, mean)
gp_025 = Cox.GammaProcess.est_survival(m3, grid_gp, x0, m->quantile(m,.025))
gp_975 = Cox.GammaProcess.est_survival(m3, grid_gp, x0, m->quantile(m,.975))
Cox.PCH.plotCI(grid_gp, [gp_025[:,2] gp_975[:,2]],col_area="orange",xlab="months")
Cox.PCH.plotCI(grid_gp, [gp_025[:,1] gp_975[:,1]],col_area=rgb(0,0,1,.5),add=true)
Cox.PCH.plotsurv(grid_gp, mean_S_gp, lwd=3, col_l=["blue","yellow"],
                 fg="grey",xlab="months", ylab="Survival Probability",
                 add=true,addlines=true);
R"lines(survfit(Surv(time,delta) ~ type, data = tongue),col='grey')";

###
println(R"coxph(Surv(time,delta) ~ type, data=tongue)")

R"pdf('../img/survival.pdf',w=13,h=7)"
R"par(mfrow=c(1,3))";
# Parametric
Cox.PCH.plotCI(t0,[param_025[:,2] param_975[:,2]],col_area="orange",
               ylab="Survival Probability",xlab="months")
Cox.PCH.plotCI(t0,[param_025[:,1] param_975[:,1]],col_area=rgb(0,0,1,.5),add=true)
Cox.PCH.plotsurv(t0,mean_S_weib,lwd=3,col_l=["blue","yellow"],add=true);
R"lines(survfit(Surv(time,delta) ~ type, data = tongue),col='grey',lwd=2)";
R"title(main='Parametric Cox Model',cex.main=2,col.main='grey30')"
R"legend('topright',legend=c('Aneuplod','Diploid','KM'),text.col=c('orange','blue','grey'),bty='n', cex=3)"

# PCH
Cox.PCH.plotCI(grid,[pch_025[:,2] pch_975[:,2]],col_area="orange",ylab="",xlab="")
Cox.PCH.plotCI(grid, [pch_025[:,1] pch_975[:,1]],col_area=rgb(0,0,1,.5),add=true)
Cox.PCH.plotsurv(grid, mean_S_pch, lwd=3, col_l=["blue","yellow"],add=true);
R"lines(survfit(Surv(time,delta) ~ type, data = tongue),col='grey',lwd=2)";
R"title(main='PCH Cox Model',cex.main=2,col.main='grey30')"

# GP
Cox.PCH.plotCI(grid_gp, [gp_025[:,2] gp_975[:,2]],col_area="orange",ylab="",xlab="")
Cox.PCH.plotCI(grid_gp, [gp_025[:,1] gp_975[:,1]],col_area=rgb(0,0,1,.5),add=true)
Cox.PCH.plotsurv(grid_gp, mean_S_gp, lwd=3, col_l=["blue","yellow"],add=true);
R"lines(survfit(Surv(time,delta) ~ type, data = tongue),col='grey',lwd=2)";
R"title(main='Gamma Process Cox Model',cex.main=2,col.main='grey30')"
R"par(mfrow=c(1,1))";
R"dev.off()"
