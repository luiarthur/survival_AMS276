immutable State
  β::Vector{Float64}
  λ::Float64
  α::Float64
  w::Vector{Float64}
  η::Float64
end

immutable Prior_β
  m::Vector{Float64}
  S::Matrix{Float64}
  cs::Matrix{Float64}
end

immutable Prior_λ
  a::Float64
  b::Float64
end

immutable Prior_α
  a::Float64
  b::Float64
  cs::Float64
end

immutable Prior_η
  a::Float64
  b::Float64
  cs::Float64
end

function fit(t::Vector{Float64}, X::Matrix{Float64}, v::Vector{Float64},
             group::Vector{Int},
             prior_β::Prior_β, prior_λ::Prior_λ,
             prior_α::Prior_α, prior_η::Prior_η,
             B::Int, burn::Int; printFreq::Int=0)

  rig(a::Float64, b::Float64) = 1 / rand(Gamma(a, 1/b))

  const N = length(unique(group)) # number of groups
  const group_ind = [ find(g->g==i, group) for i in 1:N ]
  const m = length.(group_ind)
  const S_inv = inv(prior_β.S)

  # precomputes
  const a_λ_new = prior_λ.a + sum(v)

  # log full-conditionals for Metropolis-step
  function lfc_β(β::Vector{Float64},λ::Float64,α::Float64)
    const Xb = X*β
    return -((β-prior_β.m)'S_inv*(β-prior_β.m)/2)[1] + sum(v.*Xb-λ*t.^α.*exp(Xb))
  end

  function lfc_α(α::Float64,λ::Float64,β::Vector{Float64})
    if α <= 0
      out = -Inf
    else
      out = prior_α.a+sum(m)*log(α) + α*(-prior_α.b+sum(v.*log(t))) -λ*sum(t.^α.*exp(X*β))
    end
    return out
  end
  
  function lfc_η(η::Float64,w::Vector{Float64})
    if η < 0
      out = -Inf
    else
      out = N*(η*log(η)-lgamma(η)) + η*sum(log(w)) - η*(prior_η.b+sum(w)) + (prior_η.a-1)*log(η)
    end
    return out
  end

  function update(s::State)
    const new_λ = rig(a_λ_new, prior_λ.b + sum(t.^s.α .* exp(X*s.β)) )
    #println("λ: ", new_λ)
    const new_w = [ rig(s.η + sum(v[group_ind[i]]), s.η) for i in 1:N ]
    #println("w: ", new_w)
    const new_β = MCMC.metropolis(s.β,prior_β.cs,β->lfc_β(β,new_λ,s.α))
    #println("β: ", new_β)
    const new_α = MCMC.metropolis(s.α,prior_α.cs,α->lfc_α(α,new_λ,new_β))
    #println("α: ", new_α)
    const new_η = MCMC.metropolis(s.η,prior_η.cs,η->lfc_η(η,new_w))
    #println("η: ", new_η)
    
    return State(new_β, new_λ, new_α, new_w, new_η)
  end

  const init = State(prior_β.m, #[1.93,.00832]
                     prior_λ.a/prior_λ.b,
                     prior_α.a/prior_α.b,
                     ones(N)*prior_η.a/prior_η.b, # prior for w vector
                     prior_η.a/prior_η.b)
  return MCMC.gibbs(init,update,B,burn,printFreq=printFreq)
end
