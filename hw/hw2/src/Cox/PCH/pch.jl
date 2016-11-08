immutable Priorβ
  # Prior: Normal(m,S) 
  # Proposal Step: Σ
  m::Vector{Float64}
  S::Matrix{Float64}
  Σ::Matrix{Float64}
end

immutable Priorλ
  # Prior: λⱼ~ Gamma(aⱼ, bⱼ) 
  # Proposal Step: Σ
  a::Vector{Float64}
  b::Vector{Float64}
  Σ::Matrix{Float64}
end

immutable State
  β::Vector{Float64}
  λ::Vector{Float64}
end

function pch(t::Vector{Float64}, X::Matrix{Float64}, v::Vector{Float64}, 
             grid::Vector{Float64}, priorβ::Priorβ, priorλ::Priorλ,
             B::Int, burn::Int; printFreq::Int=0)

  assert(grid[1] == 0)
  const J = length(grid) - 1
  function s(j::Int)
    assert(0 <= j <= J)
    return grid[j+1]
  end
  const (N,P) = size(X)
  assert(length(t) == N && length(v) == N)
  const S⁻¹ᵦ = inv(priorβ.S)

  function h₀(tᵢ::Float64, λ::Vector{Float64})
    assert(length(λ) == J)
    return sum([ s(j-1)<tᵢ<=s(j) ? λ[j] : 0. for j in 1:J])
  end

  function H₀(tᵢ::Float64, λ::Vector{Float64})
    assert(length(λ) == J)
    out = 0.
    j = 1
    while s(j) < tᵢ
      out += λ[j] * (s(j) - s(j-1))
      j += 1
    end
    out += λ[j] * (tᵢ - s(j-1))
    return out
  end

  function loglike(β::Vector{Float64}, λ::Vector{Float64})
    const Xb = X*β
    return sum([ v[i] * (log(h₀(t[i],λ)) + Xb[i]) - 
                 H₀(t[i],λ) * exp(Xb[i]) for i in 1:N ])
  end
  
  logpriorβ(β::Vector{Float64}) = (-(β-priorβ.m)'S⁻¹ᵦ*(β-priorβ.m)/2)[1]
  logprior_logλ(logλ::Vector{Float64}) = sum(-logλ.*(priorλ.a-1) - priorλ.b .* exp(logλ))

  function ll_plus_lp(β::Vector{Float64},logλ::Vector{Float64})
    return logpriorβ(β) + logprior_logλ(logλ) + loglike(β,exp(logλ))
  end

  function update(curr::State)
    const curr_logλ = log(curr.λ)
    const newβ = MCMC.metropolis(curr.β, priorβ.Σ, β->ll_plus_lp(β,curr_logλ))
    const new_logλ = MCMC.metropolis(curr_logλ, priorλ.Σ, logλ->ll_plus_lp(newβ,logλ))
    return State(newβ, exp(new_logλ))
  end

  const init = State(priorβ.m, priorλ.a ./ priorλ.b)
  return MCMC.gibbs(init, update, B, burn, printFreq=printFreq)
end
