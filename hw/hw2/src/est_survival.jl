function est_survival(pch::Vector{PCH.State},grid::Vector{Float64},age::Float64) 

  assert(grid[1] == 0)
  const J = length(grid) - 1
  s(j::Int) = grid[j+1]

  const B = length(pch)
  const lambda = hcat(map(p->p.λ, pch)...)'
  const beta = hcat(map(p->p.β, pch)...)'

  const P = size(beta,2)

  #const mean_lambda = vec(mean(lambda,1))
  #const mean_beta = [0.; vec(mean(beta,1))]

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
  const S = [S0[b,j] ^ exp((p==0?0:beta[b,p]) + age*beta[b,P]) for p in 0:P-1, j in 1:length(grid), b in 1:B]
  #const S = hcat([S0 .^ (exp(mean_beta[p] + age*mean_beta[P])) for p in 1:P-1]...)
  #const S = cat(3, [S0 .^ (exp(beta[:,p] + age*beta[:,P])) for p in 1:P-1]...)
  
  return S
end
