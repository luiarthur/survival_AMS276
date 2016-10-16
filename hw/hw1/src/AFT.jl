module AFT

using Distributions
import MCMC

export State_aft, aft, dic

immutable State_aft
  sig::Real
  beta::Array{Real,1}
  sig_acc::Int
  beta_acc::Array{Int,1}
end

function aft(t, X, v, init::State_aft, 
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

  logprior_beta(bj::Real,mj::Real, s2j::Real) = -(bj-mj)^2 / (2*s2j)
  logprior_sig(sig::Real) = (-a-1)*log(sig) - b/sig

  function loglike(sig::Real,beta::Array{Real,1};bj=0,jj=0)
    new_beta = collect(beta)

    if jj > 0 
      new_beta[jj] = bj
    end

    # NEED TO DO A FORMAL SPEED TEST
    #Xb = [sum(X1[i,:] .* beta) for i in 1:N] # 51s
    #Xb = sum([X1[i,j] * new_beta[j] for i in 1:N, j in 1:J],2) # 10s
    if J < 20
      if J == 2
        Xb = new_beta[1] .+ X * new_beta[2] # 2s
      else # 2 < J < 20
        Xb = sum([X1[:,j] * new_beta[j] for j in 1:J]) # 4s
      end
    else # J ≥ 20
      Xb = X1*new_beta # 121s
    end

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

  function loglike(sig::Real,beta::Array{Real,1})
    if J < 20
      if J == 2
        Xb = beta[1] .+ X * beta[2] # 2s
      else # 2 < J < 20
        Xb = sum([X1[:,j] * beta[j] for j in 1:J]) # 4s
      end
    else # J ≥ 20
      Xb = X1*beta # 121s
    end

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

  const D = [-2 * loglike(t.sig,t.beta) for t in post]
  return mean(D) + var(D) / 2
end # dic

end # AFT

# look at getfield, fieldnames
