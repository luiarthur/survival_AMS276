include("GammaProcess.jl")
include("../PCH/PCH.jl")

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

### GP
@time m3 = GammaProcess.gp(t,X,v,2.,.1,2000,1000,c=.01,printFreq=500);
R"coxph(Surv(time,delta) ~ as.factor(stage) + age, data=larynx)"
s3 = GammaProcess.summary(m3)
GammaProcess.plot(m3,"beta",[1,2,3,4]);
GammaProcess.plot(m3,"h",[1,2,3,4,5]);
GammaProcess.plot(m3,"h",[6,7,8,9,10]);

grid = sort(unique([0;t]))
age = 60.
x0 = [[0,0,0,age],[1,0,0,age],[0,1,0,age],[0,0,1,age]]
S3 = GammaProcess.est_survival(m3, grid, x0, mean)
q_lo = GammaProcess.est_survival(m3, grid, x0, m->quantile(m,.025))
q_hi = GammaProcess.est_survival(m3, grid, x0, m->quantile(m,.975))
PCH.plotsurv(grid, S3, lwd=3, col_l=["black","blue","red","green"],
             fg="grey",xlab="time", ylab="Survival Probability",
             main="Probability of Survival for different Stages",
             addlines=true);
R"lines(survfit(Surv(time,delta) ~ as.factor(stage), data=larynx))";
#rgb = R"rgb"
#PCH.plotCI(grid,[q_lo[:,1] q_hi[:,1]],col=rgb(0,0,0,.5))
#PCH.plotCI(grid,[q_lo[:,2] q_hi[:,2]],col=rgb(0,0,1,.5),add=true)
#PCH.plotCI(grid,[q_lo[:,3] q_hi[:,3]],col=rgb(1,0,0,.5),add=true)
#PCH.plotCI(grid,[q_lo[:,4] q_hi[:,4]],col=rgb(0,1,0,.5),add=true)
