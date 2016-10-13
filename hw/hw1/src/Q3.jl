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

immutable State_weib
  sig::Real
  b0::Real
  b1::Real
  sig_acc::Int
  b0_acc::Int
  b1_acc::Int
end

function update_weib(state::State_weib;m=0,s2=1,a=2,b=3,css=1,csb0=1,csb1=1)

  const sum_y = sum(y_all)

  # log-likelihood
  function loglike(sig::Real,b0::Real,b1::Real)
    sum(
        v_all .* (-log(sig) .+ log(y_all)./sig .+ (b0+b1*x_all)/sig) .- 
        exp((b0+b1)/sig) .* y_all.^(1/sig)
       )
  end

  # log-prior: beta
  logprior_beta(beta::Real) = -(beta-m)^2 / (2*s2)

  # log-prior: sig
  logprior_sig(sig::Real) = (-a=1)*log(sig) - b/sig

  # update sig
  const (new_sig, sig_acc) =
  mh_normal(state.sig, sig->loglike(sig,state.b0,state.b1)+logprior_sig(sig), 
            state.sig_acc, css, inbounds=x->x>0)

  # update b0
  const (new_b0, b0_acc) =
  mh_normal(state.b0, b0->loglike(new_sig,b0,state.b1)+logprior_beta(b0), 
            state.b0_acc, csb0)

  # update b1
  const (new_b1, b1_acc) =
  mh_normal(state.b1, b1->loglike(new_sig,new_b0,b1)+logprior_beta(b1), 
            state.b1_acc, csb1)

  return State_weib(new_sig,new_b0,new_b1,sig_acc,b0_acc,b1_acc)
end

@time out_weib = gibbs(State_weib(1,0,0,0,0,0),update_weib,B,burn);

#=
include("Q3.jl")
=#

