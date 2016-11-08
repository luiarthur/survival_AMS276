function est_survival(gp::Vector{State}, grid:: Vector{Float64}, 
                      xs::Vector{Vector{Float64}}, f)
  const S = est_survival(gp, grid, xs)
  return mapslices(s -> f(s),S,2)[:,1,:]
end

function est_survival(gp::Vector{State},grid::Vector{Float64},
                         xs::Vector{Vector{Float64}})

  return cat(3,[est_survival(gp, grid, x0) for x0 in xs]...)
end

function est_survival(gp::Vector{State},grid::Vector{Float64},x::Vector{Float64})

  assert(grid[1] == 0)
  const J = length(grid) - 1

  const B = length(gp)
  const h = map(p->p.h, gp)
  const beta = map(p->p.β, gp)

  H₀(j::Int, h::Vector{Float64}) = sum(h[1:j])

  const S = [exp(-H₀(j,h[b]) * (exp(x'(beta[b]))[1])) for j in 0:J, b in 1:B]
  
  return S
end
