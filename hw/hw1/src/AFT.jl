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
  sig_acc::Int
  beta_acc::Array{Int,1}
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
    new_beta = collect(beta)

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
    const (new_sig, sig_acc) =
      MCMC.mh_normal(state.sig, sig->loglike(sig,state.beta)+logprior_sig(sig), 
                     state.sig_acc, css, inbounds=x->x>0)
    # update bj
    const new_beta_acc = Array{Int,1}(J)
    const new_beta = collect(state.beta)
    for j in 1:J
      (new_bj, new_bj_acc) =
      MCMC.mh_normal(new_beta[j], 
                     bj -> loglike(new_sig,new_beta,bj=bj,jj=j) +
                           logprior_beta(bj,m[j],s2[j]), 
                     state.beta_acc[j], csb[j])

      new_beta[j] = new_bj
      new_beta_acc[j] = new_bj_acc
    end
    
    return State_aft(new_sig, new_beta, sig_acc, new_beta_acc)
  end # update_aft


  out = MCMC.gibbs(init, update_aft, B, burn, printFreq=printFreq)

  println()
  println("σ acceptance rate: ", out[end].sig_acc / B)
  println("β acceptance rate: ", out[end].beta_acc / B)

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
