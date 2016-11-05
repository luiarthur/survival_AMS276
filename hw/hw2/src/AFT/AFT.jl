module AFT
#=
To do:
  - Clean up code. 
  - Define all types in functino arguments
=#

using Distributions
import MCMC

include("SurvivalDistributions.jl")

export State_aft, aft, dic

immutable State_aft
  sig::Float64
  beta::Vector{Float64}
end

function loglike(sig::Float64, beta::Vector{Float64}, 
                 t::Vector{Float64}, X::Matrix{Float64}, v::Vector{Float64},
                 model::String)

  const N = length(t)
  const timeDist = let 
    if model == "lognormal"
      SurvivalDistributions.lognormal
    elseif model == "loglogistic"
      SurvivalDistributions.loglogistic
    else
      SurvivalDistributions.weibull
    end
  end

  Xb = X*beta 

  return sum([logpdf(timeDist, t[i], Xb[i], sig)*v[i] + 
              logccdf(timeDist, t[i], Xb[i], sig) *(1-v[i]) for i in 1:N])
end

function aft(t::Vector{Float64}, X::Matrix{Float64}, v::Vector{Float64},
             csb_vec::Vector{Float64},css::Float64=1.0;
             model::String="weibull",
             B::Int=1000,burn::Int=100,printFreq::Int=0)
  const P = size(X,2) + 1
  const m = zeros(Float64,P)
  const s2 = fill(1E2,P)
  const a = 2.0
  const b = 1.0
  const prior_mean_beta = m
  const prior_mean_sig = b

  aft(t,X,v,m,s2,csb_vec,a,b,css,model=model,B=B,burn=burn,printFreq=printFreq)
end

function aft(t::Vector{Float64}, X::Matrix{Float64}, v::Vector{Float64}, 
             m::Vector{Float64}, s2::Vector{Float64}, csb::Vector{Float64},
             a::Float64, b::Float64, css::Float64; # priors for σ
             model::String="weibull",
             B::Int=1000, burn::Int=100, printFreq::Int=0)

  assert(in(model,["weibull","lognormal","loglogistic"]))

  const init = State_aft(a>1 ? b/(a-1) : 1, m)
  const N = length(t)
  const X1 = [ones(N) X]
  const J = size(X1,2)
  const Σ⁻¹ᵦ = inv(Matrix(Diagonal(s2)))
  const Σ = Matrix(Diagonal(csb))

  #logprior_beta(bj::Float64,mj::Float64, s2j::Float64) = -(bj-mj)^2 / (2*s2j)
  logprior_beta(b::Vector{Float64}) = ((b-m)'Σ⁻¹ᵦ*(b-m))[1] / -2
  logprior_sig(sig::Float64) = (-a-1)*log(sig) - b/sig
  function ll(sig::Float64,beta::Vector{Float64},model::String)
    loglike(sig,beta,t,X1,v,model)
  end

  function update_aft(state::State_aft)
    # update sig
    const new_sig =
      MCMC.metropolis(state.sig, sig->ll(sig,state.beta,model)+logprior_sig(sig), 
                      css, inbounds=x->x>0)

    # Update beta
    const cand = rand(MvNormal(state.beta, Σ))
    if ll(new_sig,cand,model) + logprior_beta(cand) - 
       ll(new_sig,state.beta,model) - logprior_beta(state.beta) > log(rand())
      new_beta = cand
    else
      new_beta = copy(state.beta)
    end
    
    return State_aft(new_sig, new_beta)
  end # update_aft


  out = MCMC.gibbs(init, update_aft, B, burn, printFreq=printFreq)

  println()
  sig_acc = length(unique(map(o -> o.sig, out))) / B
  beta_acc = length(unique(map(o -> o.beta, out))) / B
  println("β acceptance rate: ", beta_acc)
  println("σ acceptance rate: ", sig_acc)

  return out
end # aft

function dic(post::Vector{State_aft}, t, X, v; model="weibull")
  assert(in(model,["weibull","lognormal","loglogistic"]))
  assert(size(X,1)==length(t))

  const N = length(t)
  const X1 = [ones(N) X]
  const J = size(X1,2)
  function ll(sig::Float64,beta::Vector{Float64},model::String)
    loglike(sig,beta,t,X1,v,model)
  end

  const D = [-2 * ll(p.sig, p.beta, model) for p in post]
  return mean(D) + var(D) / 2
end # dic

end # AFT
