function est_survival(pch::Vector{State}, grid:: Vector{Float64}, 
                      xs::Vector{Vector{Float64}}, f)
  const S = est_survival(pch, grid, xs)
  return mapslices(s -> f(s),S,2)[:,1,:]
end

function est_survival(pch::Vector{State},grid::Vector{Float64},
                         xs::Vector{Vector{Float64}})

  return cat(3,[PCH.est_survival(pch, grid, x0) for x0 in xs]...)
end

function est_survival(pch::Vector{State},grid::Vector{Float64},x::Vector{Float64})

  assert(grid[1] == 0)
  const J = length(grid) - 1
  s(j::Int) = grid[j+1]

  const B = length(pch)
  const lambda = hcat(map(p->p.λ, pch)...)'
  const beta = hcat(map(p->p.β, pch)...)'

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

  const H0 = [H₀(g, lambda[b,:]) for b in 1:B, g in grid]
  const S0 = exp(-H0)
  const S = [S0[b,j] ^ (exp(x'beta[b,:]))[1] for j in 1:length(grid), b in 1:B]
  
  return S
end
