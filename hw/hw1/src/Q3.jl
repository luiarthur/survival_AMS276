using RCall, MCMC, Distributions
srand(256)

R"""
library(fields)
library(KMsurv)   # To get the datasets in K-M
library(survival) # R functions
library(rcommon)  # devtools::install_github('luiarthur/rcommon')
data(tongue)
tongue$time <- tongue$time/10 # to obtain computational stability
""";

#R"my.pairs(tongue)";
#R"my.pairs(cbind(tongue[,1],log(tongue[,2]),tongue[,3]))";

R"tongue$type <- ifelse(tongue$type==1,1,0)"
R"aneuploid <- tongue[which(tongue$type==1),]" 
R"diploid <- tongue[which(tongue$type!=1),]" 

@rget aneuploid diploid tongue

###


# Precomputes:
const y_a = aneuploid[:,2]
const v_a = aneuploid[:,3]

const y_d = diploid[:,2]
const v_d = diploid[:,3]

const y_all = tongue[:,2]
const v_all = tongue[:,3]
const x_all = tongue[:,1]

immutable State
  lambda::Real
  alpha::Real
end
                           # Priors:
function update_gen(y,v,cs;a_alpha=.1,b_alpha=.1,a_lam=0,b_lam=0)
  const sum_v = sum(v)

  function update(state::State)
    #update lambda
    const new_lambda = rand(Gamma( a_lam+sum_v, 1/(sum(y.^state.alpha)+b_lam) ))

    #update alpha (mh)
    function loglike_plus_logprior(alpha::Real)
      (a_alpha+sum_v-1) * log(alpha) + 
      sum(alpha.*v.*log(y)-new_lambda.*y.^alpha) -alpha*b_alpha
    end

    const new_alpha = 
      metropolis(state.alpha, loglike_plus_logprior, cs, inbounds=x->x>0)

      State(new_lambda,new_alpha)
  end

  return update
end

const B = 10000
const burn = 1000

update_a = update_gen(y_a,v_a,.3;a_alpha=.1,b_alpha=.1,a_lam=0,b_lam=0)
update_d = update_gen(y_d,v_d,.3;a_alpha=.1,b_alpha=.1,a_lam=0,b_lam=0)
@time out_a = gibbs(State(1,1),update_a,B,burn);
@time out_d = gibbs(State(1,1),update_d,B,burn);

lam_a= map(o->o.lambda, out_a)
alp_a = map(o->o.alpha, out_a)
acc_a = length(unique(alp_a)) / length(alp_a)

lam_d= map(o->o.lambda, out_d)
alp_d = map(o->o.alpha, out_d)
acc_d = length(unique(alp_d)) / length(alp_d)

@rput lam_a alp_a lam_d alp_d;
R"pdf('../img/post_a.pdf')"
R"plotPosts(cbind(lam_a,alp_a),c('lambda','alpha'),legend.pos='right',cex.l=1.3,show.x=F,cex.a=2)";
R"dev.off()"
R"pdf('../img/post_b.pdf')"
R"plotPosts(cbind(lam_d,alp_d),c('lambda','alpha'),legend.pos='right',cex.l=1.3,show.x=F,cex.a=2)";
R"dev.off()"

# Compare
surv_weib(y,a,b) = exp(-b*y^a)
S_a(y) = map(o -> surv_weib(y,o.alpha,o.lambda), out_a)
S_d(y) = map(o -> surv_weib(y,o.alpha,o.lambda), out_d)
yseq = linspace(0,50,100)
Sa = map(o -> S_a(o), yseq)
Sd = map(o -> S_d(o), yseq)

@rput Sa Sd yseq

R"pdf('../img/and.pdf')"
R"plot(0,xlim=range(yseq),ylim=0:1,cex=0,main='Survival Probability of Aneuploid and Diploid Tumor Patients',ylab='Survival Probability',xlab='time (10 weeks)',bty='n',fg='grey')"
R"color.btwn(yseq,sapply(Sa,quantile,.025),sapply(Sa,quantile,.975),from=0,to=50,col=rgb(0,0,1,.3))"
R"color.btwn(yseq,sapply(Sd,quantile,.025),sapply(Sd,quantile,.975),from=0,to=50,col=rgb(1,0,0,.3))"
R"legend('topright',legend=c('Aneuploid','Diploid'),text.col=c('dodgerblue','pink'),text.font=2,bty='n',cex=3)"
R"dev.off()"

# Compare alpha
R"pdf('../img/compareAlpha.pdf')"
R"color.den(density(alp_a),from=0,to=1.3,col.area=rgb(0,0,1,.4),add=F,fg='grey',bty='n',main='')"
R"color.den(density(alp_d),from=0,to=1.3,col.area=rgb(1,0,0,.4),add=T,col.den=rgb(0,0,0,0))"
R"legend('topright',legend=c('Aneuploid','Diploid'),text.col=c('dodgerblue','pink'),text.font=2,bty='n',cex=2.6)"
R"dev.off()"

R"pdf('../img/alphadiff.pdf')"
R"plotPost(alp_a-alp_d,main=expression(alpha[diff]),cex.main=3,cex.l=2,cex.a=2)"
R"dev.off()"

# Sensitivity on λ and α
function sensitivity(a_alpha,b_alpha,a_lam,b_lam; cs=.3)
  update = update_gen(y_a,v_a,cs;a_alpha=a_alpha,b_alpha=b_alpha,a_lam=a_lam,b_lam=b_lam)
  out_sens = gibbs(State(1,1,0),update,500,1000);
  lam_sens = map(o -> o.lambda, out_sens)
  alp_sens = map(o->o.alpha, out_sens)
  Dict(
       :a=>Dict(:mean=>mean(alp_sens),:ci=>quantile(alp_sens,[.025,.975])),
       :l=>Dict(:mean=>mean(lam_sens),:ci=>quantile(lam_sens,[.025,.975])),
      )
end

# for lambda improper
n = 1000
U = rand(Uniform(0,5),(n,2))
@time tmp=sensitivity(.1,.1,0,0);
sens_a = Array{typeof(tmp),1}(n)
sens_l = Array{typeof(tmp),1}(n)

# 90s on 4core, o.w. double on 1 thread
@time Threads.@threads for i in 1:n
  sens_a[i] = sensitivity(U[i,1],U[i,2],0,0)
  sens_l[i] = sensitivity(.1,.1,U[i,1],U[i,2])
end

function plotSens(sens,param="a")
  a_means = map(o -> o[:a][:mean],sens)
  l_means = map(o -> o[:l][:mean],sens)
  @rput U a_means l_means param
  R"""
  par(mfrow=c(1,2))
  par(mar=c(4,8,1,5))
  xlab=ifelse(param=="a",expression(alpha[shape]),expression(lambda[shape]))
  ylab=ifelse(param=="a",expression(alpha[rate]),expression(lambda[rate]))
  quilt.plot(U[,1],U[,2],a_means,xlab=xlab,ylab=ylab,cex.main=2,
             fg='grey',bty='n',cex.lab=2,
             main=expression(alpha))
  quilt.plot(U[,1],U[,2],l_means,xlab=xlab,cex.main=2,
             fg='grey',bty='n',cex.lab=2,
             main=expression(lambda))
  par(mfrow=c(1,1),las=0)
  """
end
R"pdf('../img/sensa.pdf',w=10,h=5)"
plotSens(sens_a,"a")
R"dev.off()"
R"pdf('../img/sensl.pdf',w=10,h=5)"
plotSens(sens_l,"l")
R"dev.off()"

#################################
srand(256)
#AFT Models
include("AFT.jl")
const N = length(y_all)
y = collect(y_all)
X = collect(x_all')'
v = collect(v_all * 1.0)

# Weibull
@time aft_weib = AFT.aft(y, X, v, B=B, burn=1000);

aft_weib_sig = map(o -> o.sig, aft_weib)
aft_weib_b0 = map(o -> o.beta[1], aft_weib)
aft_weib_b1 = map(o -> o.beta[2], aft_weib)

@rput aft_weib_sig aft_weib_b0 aft_weib_b1
R"pdf('../img/aft_weib.pdf')"
R"plotPosts(cbind(aft_weib_sig, aft_weib_b0, aft_weib_b1),legend.pos='right',cex.l=1,show.x=F,cex.a=1.5)"; println()
R"dev.off()"

# LogLogistic
@time aft_loglog = AFT.aft(y, X, v, B=B, burn=1000, model="loglogistic");

aft_loglog_sig = map(o -> o.sig, aft_loglog)
aft_loglog_b0 = map(o -> o.beta[1], aft_loglog)
aft_loglog_b1 = map(o -> o.beta[2], aft_loglog)

@rput aft_loglog_sig aft_loglog_b0 aft_loglog_b1
R"pdf('../img/aft_loglog.pdf')"
R"plotPosts(cbind(aft_loglog_sig, aft_loglog_b0, aft_loglog_b1),legend.pos='right',cex.l=1,show.x=F,cex.a=1.5)"; println()
R"dev.off()"

# LogNormal
@time aft_logNorm = AFT.aft(y, X, v, B=B, burn=1000, model="lognormal");

aft_logNorm_sig = map(o -> o.sig, aft_logNorm)
aft_logNorm_b0 = map(o -> o.beta[1], aft_logNorm)
aft_logNorm_b1 = map(o -> o.beta[2], aft_logNorm)

@rput aft_logNorm_sig aft_logNorm_b0 aft_logNorm_b1
R"pdf('../img/aft_lognorm.pdf')"
R"plotPosts(cbind(aft_logNorm_sig, aft_logNorm_b0, aft_logNorm_b1),legend.pos='right',cex.l=1,show.x=F,cex.a=1.5)"; println()
R"dev.off()"

### Print R models:
println("Equivalent models in R")
R"weib_mod <- survreg(Surv(time, delta) ~ as.factor(type), dist='weibull', data=tongue)"
R"loglogistic_mod <- survreg(Surv(time, delta) ~ as.factor(type), dist='loglogistic', data=tongue)"
R"lognormal_mod <- survreg(Surv(time, delta) ~ as.factor(type), dist='lognormal', data=tongue)" 

R"print(summary(weib_mod))"
R"print(summary(loglogistic_mod))"
R"print(summary(lognormal_mod))"

### DIC:
weib_dic =    AFT.dic(aft_weib,y,X,v, model="weibull")
loglog_dic =  AFT.dic(aft_loglog,y,X,v, model="loglogistic")
lognorm_dic = AFT.dic(aft_logNorm,y,X,v, model="lognormal")

println("DIC for Weibull: ",weib_dic)
println("DIC for loglog:  ",loglog_dic)
println("DIC for lognorm: ",lognorm_dic)

### Posterior Acceleration (to death) Factor:
println()
weibAF = map(o -> exp(o.beta[2]), aft_weib)
println("Weibull Acceleration Factor (for aneploid compared to diploid): ",mean(weibAF))
@rput weibAF
R"pdf('../img/weibaf.pdf')"
R"plotPost(weibAF,main='',legend.pos='right',cex.l=1.5)"
R"dev.off()"

loglogAF = map(o -> exp(o.beta[2]), aft_loglog)
println("Loglog Acceleration Factor (for aneploid compared to diploid): ",mean(loglogAF))
@rput loglogAF
R"pdf('../img/loglogaf.pdf')"
R"plotPost(loglogAF,main='',legend.pos='right',cex.l=1.5)"
R"dev.off()"

lognormAF = map(o -> exp(o.beta[2]), aft_logNorm)
println("LogNorm Acceleration Factor (for aneploid compared to diploid): ",mean(lognormAF))
@rput lognormAF
R"pdf('../img/lognormaf.pdf')"
R"plotPost(lognormAF,main='',legend.pos='right',cex.l=1.5)"
R"dev.off()"

#=
include("Q3.jl")
=#
