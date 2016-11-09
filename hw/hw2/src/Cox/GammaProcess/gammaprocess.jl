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

  # Can be precomputed
  function part(j::Int) # (safe, death) sets
    const risk = Vector{Int}() # Risk set
    #const safe  = Vector{Int}()  # Safe set: At risk but not failed (Rⱼ - Dⱼ)
    const death  = Vector{Int}()  # Death set: At risk and failed (Dⱼ)
    for i in 1:N
      if t[i] > s(j-1)
        push!(risk,i)
      end
      if s(j-1) < t[i] <= s(j) && v[i] == 1
        push!(death,i)
      end
    end
    const safe = collect(setdiff(risk,death))
    return (safe,death)
  end

  const safe = [part(j)[1] for j in 1:J]
  const death = [part(j)[2] for j in 1:J]
  #return (safe,death)

  function loglike(β::Vector{Float64}, h::Vector{Float64}, j::Int)
    const Xb = X*β
    #const (safe, death) = part(j)
    -h[j] * sum([exp(Xb[i]) for i in safe[j]]) + # log survival
    sum([log( 1-exp(-h[j]*exp(Xb[i])) ) for i in death[j]]) # log hazard
  end

  function loglike(β::Vector{Float64}, h::Vector{Float64})
    return sum([loglike(β, h, j) for j in 1:J])
  end
  
  logpriorβ(β::Vector{Float64}) = (-(β-priorᵦ.m)'S⁻¹ᵦ*(β-priorᵦ.m)/2)[1]
  function logprior_logh(logh::Vector{Float64}) 
    return sum(-logh.*(priorₕ.ac-1) - priorₕ.c * exp(logh))
  end

  function ll_plus_lp(β::Vector{Float64},logh::Vector{Float64})
    return logpriorβ(β) + logprior_logh(logh) + loglike(β,exp(logh))
  end

  function update(curr::State)
    const curr_logh = log(curr.h)
    const newβ = MCMC.metropolis(curr.β, priorᵦ.Σ, β->ll_plus_lp(β,curr_logh))
    const new_logh = MCMC.metropolis(curr_logh, priorₕ.Σ, logh->ll_plus_lp(newβ,logh))
    return State(newβ, exp(new_logh))
  end

  const init = State(priorᵦ.m, priorₕ.ac / priorₕ.c)
  return MCMC.gibbs(init, update, B, burn, printFreq=printFreq)
end # gp

function gp(t::Vector{Float64}, X::Matrix{Float64}, v::Vector{Float64}, 
            csᵦ::Float64, csₕ::Float64,
            B::Int, burn::Int; η::Float64=1., c::Float64=1., κ::Float64=1., printFreq::Int=0)
  const grid = sort(unique([0;t]))
  const J = length(grid)-1
  const a = (grid[2:end].^κ - grid[1:J].^κ) * η
  const (N,P) = size(X)
  const Σᵦ = sym(inv(X'X))
  const priorᵦ = Priorᵦ(fill(0.,P), eye(P)*100., Σᵦ*csᵦ)
  const priorₕ = Priorₕ(c, a*c, eye(J) * csₕ)

  return gp(t,X,v,grid, priorᵦ, priorₕ, B, burn, printFreq=printFreq)
end
