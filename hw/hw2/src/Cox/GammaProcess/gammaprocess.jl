sym(M::Matrix{Float64}) = (M' + M) / 2
immutable Priorᵦ
  # Prior: Normal(m,S) 
  # Proposal Step: Σ
  m::Vector{Float64}
  S::Matrix{Float64}
  Σ::Matrix{Float64}
end

immutable Priorₕ
  # Prior: hⱼ~ Gamma(c₀[H*(sⱼ) - H*(sⱼ-₁)], c₀)
  # Proposal Step: Σ
  c::Float64
  ac::Vector{Float64} # J intervals
  Σ::Matrix{Float64}
end

immutable State
  β::Vector{Float64}
  h::Vector{Float64}
end

function gp(t::Vector{Float64}, X::Matrix{Float64}, v::Vector{Float64}, 
            grid::Vector{Float64}, priorᵦ::Priorᵦ, priorₕ::Priorₕ,
            B::Int, burn::Int; printFreq::Int=0)

  assert(grid[1] == 0 && grid[end] >= maximum(t))
  const (N,P) = size(X)
  assert(length(t) == N && length(v) == N)
  const J = length(grid) - 1

  s(j::Int) = grid[j+1] # shift idx

  const S⁻¹ᵦ = inv(priorᵦ.S)

  function part(j::Int) # (safe, death) sets
    const risk = Vector{Int}()   # Risk set: At risk as of interval j
    const death  = Vector{Int}() # Death set: At risk and failed (Dⱼ)
    for i in 1:N
      if t[i] > s(j-1)
        push!(risk,i)
      end
      if s(j-1) < t[i] <= s(j) && v[i] == 1
        push!(death,i)
      end
    end
    const safe = collect(setdiff(risk,death)) # Safe set: At risk but not failed (Rⱼ - Dⱼ)
    return (safe,death)
  end

  const safe = [part(j)[1] for j in 1:J]
  const death = [part(j)[2] for j in 1:J]

  function loglike(β::Vector{Float64}, h::Vector{Float64})
    const Xb = X*β
    out = 0.
    for j in 1:J
      out += -h[j]*sum(exp(Xb[safe[j]])) + sum(log(1-exp(-h[j]*exp(Xb[death[j]]))))
    end
    return out
  end
  
  logpriorβ(β::Vector{Float64}) = (-(β-priorᵦ.m)'S⁻¹ᵦ*(β-priorᵦ.m)/2)[1]
  function logprior_logh(logh::Vector{Float64}) 
    return sum(logh.*priorₕ.ac - priorₕ.c*exp(logh))
  end

  function update(curr::State)
    const curr_logh = log(curr.h)
    const newβ = MCMC.metropolis(curr.β, priorᵦ.Σ, β-> loglike(β,curr.h) + logpriorβ(β))
    const new_logh = MCMC.metropolis(curr_logh, priorₕ.Σ, logh->loglike(newβ,exp(logh)) + logprior_logh(logh))
    return State(newβ, exp(new_logh))
  end

  const init = State(priorᵦ.m, priorₕ.ac / priorₕ.c)

  # Diagnostics: Probably should :ard code metropolis update
  #@time loglike(init.β,init.h)
  #@time loglike(init.β,init.h)
  #@time logpriorβ(init.β)
  #@time logpriorβ(init.β)
  #@time logprior_logh(log(init.h))
  #@time logprior_logh(log(init.h))
  return MCMC.gibbs(init, update, B, burn, printFreq=printFreq)
end # gp

"""
Gamma process

    gp(t::Vector{Float64}, X::Matrix{Float64}, v::Vector{Float64}, 
                csβ::Float64, csh::Float64,
                B::Int, burn::Int; 
                η::Float64=1., c::Float64=1., κ::Float64=1., printFreq::Int=0)


# Arguments:
* `t`:     time vector
* `X`:     covariate matrix
* `v`:     indicator for observed failure (1) or right-censorship (0)
* `grid`:  grid locations to do nonparametric estimates. For example, `grid = sort(unique([0; t]))`, which is the default.
* `csβ`:   proposal step-size (scale) for β vector (increase if accβ is too high, decrease if too low)
* `csh`:   proposal step-size (scale) for h vector (increase if acch is too high, decrease if too low)
* `B`:     number of mcmc samples
* `burn`:  burn-in
* `c,η,κ`: The baseline cumulative hazard H₀(t)∼ GammaProcess(c(ηtᵏ), c). 
  * smaller `c` reflects higher prior uncertainty around the baseline cumulative hazard, which is Weibull.
  * `c,η,κ` each have positive support
"""
function gp(t::Vector{Float64}, X::Matrix{Float64}, v::Vector{Float64}, 
            grid::Vector{Float64}, csβ::Float64, csh::Float64, 
            B::Int, burn::Int; η::Float64=1., c::Float64=1., κ::Float64=1., printFreq::Int=0)

  const J = length(grid)-1
  const a = (grid[2:end].^κ - grid[1:J].^κ) * η
  const (N,P) = size(X)
  const Σᵦ = sym(inv(X'X))
  const priorᵦ = Priorᵦ(fill(0.,P), eye(P)*100., Σᵦ*csβ)
  const priorₕ = Priorₕ(c, a*c, eye(J) * csh)

  return gp(t,X,v,grid, priorᵦ, priorₕ, B, burn, printFreq=printFreq)
end

function gp(t::Vector{Float64}, X::Matrix{Float64}, v::Vector{Float64}, 
            csβ::Float64, csh::Float64,
            B::Int, burn::Int; η::Float64=1., c::Float64=1., κ::Float64=1., printFreq::Int=0)
  const grid = sort(unique([0; t]))
  gp(t,X,v,grid,csβ, csh, B, burn, η=η, c=c, κ=κ, printFreq=printFreq)
end

