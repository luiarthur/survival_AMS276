module AFT

using Distributions
import MCMC

export State_weib, aft

immutable State_weib
  sig::Real
  beta::Array{Real,1}
  sig_acc::Int
  beta_acc::Array{Int,1}
end

function aft(y, X, v, init::State_weib, 
             m, s2, csb, # priors for β
             a, b, css; # priors for σ
             B=1000, burn=100, printFreq=0)

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
    # For J = 2,
    #Xb = [sum(X1[i,:] .* beta) for i in 1:N] # 51s
    #Xb = sum([X1[i,j] * new_beta[j] for i in 1:N, j in 1:J],2) # 10s
    if J < 20
      Xb = sum([X1[:,j] * new_beta[j] for j in 1:J]) # 4s
    else 
      Xb = X1*new_beta # 121s
    end
    sum(v .* (-log(sig) .+ log(y)/sig .+ Xb/sig) .- exp(Xb/sig) .* y.^(1/sig))
  end

  function update_weib(state::State_weib)
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
    
    return State_weib(new_sig, new_beta, sig_acc, new_beta_acc)
  end # update weibull

  return MCMC.gibbs(init, update_weib, B, burn, printFreq=printFreq)
end # aft

end # AFT
