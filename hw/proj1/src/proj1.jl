println("Starting...")
include("Frailty/Frailty.jl")
using RCall
R"library(rcommon)"
plotpost = R"plotPosts"
sym(M::Matrix{Float64}) = (M' + M) / 2
srand(276);

#=
model:
# h(tᵢⱼ| wᵢ,xᵢⱼ) = h₀(tᵢⱼ) wᵢ exp(xᵢⱼ'β)
=#
R"""
library(survival)
kidney <- read.table("dat/kidney.txt",header=TRUE,skip=5)
kidney$sex <- ifelse(kidney$sex==2, 1, 0) # M=0, F=1
""";
@rget kidney;

const t = kidney[:time] * 1.;
const X = [kidney[:age] kidney[:sex]] * 1.; # age, sex
const v = kidney[:nu] * 1.;
const group = convert(Vector{Int64},kidney[:cluster]);

prior_β = Frailty.Prior_β(zeros(2), eye(2)*10^3, sym(inv(X'X))*5.)
prior_λ = Frailty.Prior_λ(.001,.001)
prior_α = Frailty.Prior_α(.001,.001,.1)
prior_η = Frailty.Prior_η(.001,.001,2)

@time out = Frailty.fit(t,X,v,group,prior_β,prior_λ,prior_α,prior_η,10000,1000);

model = Frailty.summary(out)
println(R"coxph(Surv(time,nu) ~ age+sex+frailty(cluster,distribution='gaussian'), data=kidney)")
println(model)

R"pdf('../img/beta.pdf')"
plotpost(hcat(map(m->m.β,out)...)',cnames=["age","sex"]);
R"dev.off()"

R"pdf('../img/gamma.pdf')"
plotpost(map(m->m.λ,out),main="gamma");
R"dev.off()"
R"pdf('../img/alpha.pdf')"
plotpost(map(m->m.α,out),main="alpha");
R"dev.off()"
R"pdf('../img/eta.pdf')"
plotpost(map(m->1/m.η,out),main="kappa = 1/eta");
R"dev.off()"

w_ci = model.w.q
w_mean = model.w.MEAN
@rput w_mean w_ci;
R"pdf('../img/w.pdf')"
R"""
tmp_N <- length(w_mean)
plot(w_mean,tmp_N:1,xlim=c(0,4.5),pch=20,col="dodgerblue",cex=3,
     yaxt='n',bty='n',fg='grey',xlab='Frailty',ylab='Cluster',
     main="Frailty by Cluster")
add.errbar(ci=w_ci,transpose=TRUE,x=tmp_N:1,col="dodgerblue",lwd=1,lty=2)
axis(2,at=1:tmp_N,label=tmp_N:1,las=2,cex.axis=.8,col='grey')
""";
R"dev.off()"

#=
include("proj1.jl")
=#
