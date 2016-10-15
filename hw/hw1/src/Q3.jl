using RCall, MCMC, Distributions

R"""
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
  alpha_acc::Int
end
                           # Priors:
function update_gen(y,v,cs;a_alpha=2,b_alpha=1,a_lam=0,b_lam=0)
  const sum_v = sum(v)

  function update(state::State)
    #update lambda
    const new_lambda = rand(Gamma( a_lam+sum_v, 1/(sum(y.^state.alpha)+b_lam) ))

    #update alpha (mh)
    function loglike_plus_logprior(alpha::Real)
      (a_alpha+sum_v-1) * log(alpha) + 
      sum(alpha.*v.*log(y)-new_lambda.*y.^alpha) -alpha*b_alpha
    end

    const (new_alpha,new_alpha_acc) = 
      mh_normal(state.alpha, loglike_plus_logprior, state.alpha_acc, cs, 
                inbounds=x->x>0)

    State(new_lambda,new_alpha, new_alpha_acc)
  end

  return update
end

const B = 10000
const burn = 1000

update_a = update_gen(y_a,v_a,.3;a_alpha=2,b_alpha=1,a_lam=0,b_lam=0)
update_d = update_gen(y_d,v_d,.3;a_alpha=2,b_alpha=1,a_lam=0,b_lam=0)
@time out_a = gibbs(State(1,1,0),update_a,B,burn);
@time out_d = gibbs(State(1,1,0),update_d,B,burn);

lam_a= map(o->o.lambda, out_a)
alp_a = map(o->o.alpha, out_a)
acc_a = out_a[end].alpha_acc / length(out_a)

lam_d= map(o->o.lambda, out_d)
alp_d = map(o->o.alpha, out_d)
acc_d = out_d[end].alpha_acc / length(out_d)

@rput lam_a alp_a lam_d alp_d;
R"plotPosts(cbind(lam_a,alp_a),c('lambda','alpha'),legend.pos='right')";
R"plotPosts(cbind(lam_d,alp_d),c('lambda','alpha'),legend.pos='right')";

#################################

#AFT Models
include("AFT.jl")
init = AFT.State_aft(1,zeros(2),1,zeros(Int,2))
const N = length(y_all)

# Weibull
@time aft_weib = AFT.aft(y_all, x_all, v_all, init, 
                         zeros(2), ones(2), [1,1], 
                         2,3,1,printFreq=10,B=B,burn=5000);

aft_weib_sig = map(o -> o.sig, aft_weib)
aft_weib_b0 = map(o -> o.beta[1], aft_weib)
aft_weib_b1 = map(o -> o.beta[2], aft_weib)

@rput aft_weib_sig aft_weib_b0 aft_weib_b1
R"plotPosts(cbind(aft_weib_sig, aft_weib_b0, aft_weib_b1),legend.pos='right',cex.l=.8,show.x=F)"; println()

# LogLogistic
@time aft_loglog = AFT.aft(y_all, x_all, v_all, init, 
                           zeros(2), ones(2), [1,1], 
                           2,3,1,printFreq=10,B=B,burn=5000,model="loglogistic");

aft_loglog_sig = map(o -> o.sig, aft_loglog)
aft_loglog_b0 = map(o -> o.beta[1], aft_loglog)
aft_loglog_b1 = map(o -> o.beta[2], aft_loglog)

@rput aft_loglog_sig aft_loglog_b0 aft_loglog_b1
R"plotPosts(cbind(aft_loglog_sig, aft_loglog_b0, aft_loglog_b1),legend.pos='right',cex.l=.8,show.x=F)"; println()

# LogNormal
@time aft_logNorm = AFT.aft(y_all, x_all, v_all, init, 
                            zeros(2), ones(2), [1,1], 
                            2,3,1,printFreq=10,B=B,burn=5000,model="lognormal");

aft_logNorm_sig = map(o -> o.sig, aft_logNorm)
aft_logNorm_b0 = map(o -> o.beta[1], aft_logNorm)
aft_logNorm_b1 = map(o -> o.beta[2], aft_logNorm)

@rput aft_logNorm_sig aft_logNorm_b0 aft_logNorm_b1
R"plotPosts(cbind(aft_logNorm_sig, aft_logNorm_b0, aft_logNorm_b1),legend.pos='right',cex.l=.8,show.x=F)"; println()

#=
include("Q3.jl")
R"weib_mod <- survreg(Surv(time, delta) ~ as.factor(type), dist='weibull', data=tongue)"
R"loglogistic_mod <- survreg(Surv(time, delta) ~ as.factor(type), dist='loglogistic', data=tongue)"
R"lognormal_mod <- survreg(Surv(time, delta) ~ as.factor(type), dist='lognormal', data=tongue)" 
=#
