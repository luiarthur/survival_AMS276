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
end

function fit(t::Vector{Float64}, X::Matrix{Float64}, v::Vector{Float64},
             group::Vector{Int},
             prior_β::Prior_β, prior_λ::Prior_λ,
             prior_α::Prior_α, prior_η::Prior_η, B::Int, burn::Int;
             printFreq::Int=0)

  rig(a::Float64, b::Float64) = rand(Gamma(a, 1/b))

  const N = length(unique(group)) # number of groups
  const group_ind = [ find(g->g==i, group) for i in 1:N ]
  const m = length.(group_ind)

  # precomputes
  const a_λ_new = prior_λ.a + sum(v)

  function update(s::State)
    const new_λ = rig(a_λ_new, prior_λ.b * sum(t.^s.α .* X*s.β) )
    const new_η = rig(prior_η.a, prior_η.b + sum(m.*(s.w-log(s.w))) )
    const new_w = [ rig(new_η + sum(v[group_ind[i]]), new_η) for i in 1:N ]

    const new_β = MCMC.metropolis()
    const new_α = MCMC.metropolis()
  end

  const init = State(prior_β.m,
                     prior_λ.a/prior_λ.b
                     prior_α.a/prior_α.b,
                     ones(N)*prior_η.a/prior_η.b, # prior for w vector
                     prior_η.a/prior_η.b)
  return MCMC.gibbs(init,update,B,burn,printFreq=printFreq)
end
