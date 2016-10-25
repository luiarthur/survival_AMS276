module AFT
#=
To do:
  - Clean up code. 
  - Define all types in functino arguments
=#

using Distributions
import MCMC

export State_aft, aft, dic

immutable State_aft
  sig::Float64
  beta::Array{Float64,1}
end

function aft(t::Array{Float64,1}, X::Array{Float64,2}, v::Array{Float64,1};
             csb=1.0,css=1.0,model="weibull",B=1000,burn=100,printFreq=0)
  const P = size(X,2) + 1
  const m = zeros(Float64,P)
  const s2 = fill(10.0,P)
  const a = 2.0
  const b = 1.0
  const prior_mean_beta = m
  const prior_mean_sig = b
  const init = State_aft(prior_mean_sig, prior_mean_beta)
  const csb_vec = fill(csb,P)

  aft(t,X,v,init,m,s2,csb_vec,a,b,css,model=model,B=B,burn=burn,
      printFreq=printFreq)
end

function aft(t::Array{Float64,1}, X::Array{Float64,2}, v::Array{Float64,1}, 
             init::State_aft, 
             m, s2, csb, # priors for β
             a, b, css; # priors for σ
             model="weibull",
             B=1000, burn=100, printFreq=0)

  assert(in(model,["weibull","lognormal","loglogistic"]))

  
  const y = in(model,["loglogistic","lognormal"]) ? log(t) : t
  const sum_y = sum(y)
  const N = length(y)
  const X1 = [ones(N) X]
  const J = size(X1,2)

  logprior_beta(bj::Float64,mj::Float64, s2j::Float64) = -(bj-mj)^2 / (2*s2j)
  logprior_sig(sig::Float64) = (-a-1)*log(sig) - b/sig

  function loglike(sig::Float64,beta::Array{Float64,1};bj=0,jj=0)
    new_beta = copy(beta)

    if jj > 0 
      new_beta[jj] = bj
    end

    Xb = X1*new_beta # This is a case where defining types really speeds up things

    if model == "lognormal"
      lam = (y .+ Xb) / sig
      sum(v .* -(log(sig) .+ (lam.^2)/2) + 
          (1 .- v) .* log([1-cdf(Normal(0,1),l) for l in lam]) )
    elseif model == "loglogistic"
      lam = (y .+ Xb) / sig
      sum(v .* (lam .- log(sig) - log(1.+exp(lam))) .- log(1.+exp(lam)))
    else
      sum(v .* (-log(sig) .+ log(y)/sig .+ Xb/sig) .- exp(Xb/sig) .* y.^(1/sig))
    end
  end

  function update_aft(state::State_aft)
    # update sig
    const new_sig =
      MCMC.metropolis(state.sig, sig->loglike(sig,state.beta)+logprior_sig(sig), 
                      css, inbounds=x->x>0)
    # update bj
    const new_beta = copy(state.beta)
    for j in 1:J
      new_bj =
      MCMC.metropolis(new_beta[j], 
                      bj -> loglike(new_sig,new_beta,bj=bj,jj=j) +
                            logprior_beta(bj,m[j],s2[j]), 
                      csb[j])

      new_beta[j] = new_bj
    end
    
    return State_aft(new_sig, new_beta)
  end # update_aft


  out = MCMC.gibbs(init, update_aft, B, burn, printFreq=printFreq)

  println()
  sig_acc = length(unique(map(o -> o.sig, out))) / B
  beta_acc = length(unique(map(o -> o.beta, out))) / B
  println("σ acceptance rate: ", sig_acc)
  println("β acceptance rate: ", beta_acc)

  return out
end # aft

function dic(post::Array{State_aft,1}, t, X, v; model="weibull")
  assert(in(model,["weibull","lognormal","loglogistic"]))
  assert(size(X,1)==length(t))

  const y = in(model,["loglogistic","lognormal"]) ? log(t) : t
  const N = length(y)
  const X1 = [ones(N) X]
  const J = size(X1,2)

  function loglike(sig::Float64,beta::Array{Float64,1})
    Xb = X1*beta

    if model == "lognormal"
      lam = (y .+ Xb) / sig
      sum(v .* -(.5*log(2*pi*sig^2) .+ (lam.^2)/2) + 
          (1 .- v) .* log([1-cdf(Normal(0,1),l) for l in lam]) )
    elseif model == "loglogistic"
      lam = (y .+ Xb) / sig
      sum(v .* (lam .- log(sig) - log(1.+exp(lam))) .- log(1.+exp(lam)))
    else
      sum(v .* (-log(sig) .+ log(y)/sig .+ Xb/sig) .- exp(Xb/sig) .* y.^(1/sig))
    end
  end

  const D = [-2 * loglike(t.sig,t.beta) for t in post]
  return mean(D) + var(D) / 2
end # dic

end # AFT

# look at getfield, fieldnames
