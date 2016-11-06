immutable Priorβ
  m::Vector{Float64}
  S::Matrix{Float64}
  Σ::Matrix{Float64}
end

immutable Priorλ
  m::Vector{Float64}
  S::Matrix{Float64}
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
  const S⁻¹ₗ = inv(priorλ.S)

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
    if any(λ .< 0) 
      out = 0. 
    else
      out = sum([ v[i] * (log(h₀(t[i],λ)) + Xb[i]) - 
                 (1-v[i]) * H₀(t[i],λ) * exp(Xb[i]) for i in 1:N ])
    end

    return out
  end
  
  logpriorβ(β::Vector{Float64}) = (-(β-priorβ.m)'S⁻¹ᵦ*(β-priorβ.m)/2)[1]
  function logpriorλ(λ::Vector{Float64})
    if any(λ .< 0)
      out = -Inf
    else
      out = sum((.1-1) * log(λ) - .1 * λ)
    end
    return out
  end
  #logpriorλ(λ::Vector{Float64}) = (-(λ-priorλ.m)'S⁻¹ₗ*(λ-priorλ.m)/2)[1]
  function ll_plus_lp(β::Vector{Float64},λ::Vector{Float64})
    #return logpriorβ(β) + logpriorλ(λ) + loglike(β,exp(λ))
    return logpriorβ(β) + logpriorλ(λ) + loglike(β,λ)
  end

  function update(curr::State)
    const newβ = MCMC.metropolis(curr.β, priorβ.Σ, β->ll_plus_lp(β,curr.λ))
    const newλ = MCMC.metropolis(curr.λ, priorλ.Σ, λ->ll_plus_lp(newβ,λ))
    return State(newβ, newλ)
  end

  const init = State(priorβ.m, priorλ.m)
  return MCMC.gibbs(init, update, B, burn, printFreq=printFreq)
end
