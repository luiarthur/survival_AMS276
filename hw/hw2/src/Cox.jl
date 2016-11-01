module Cox

using Distributions
import MCMC.gibbs, Base.show

export coxph_weibull, summary_cox

function part(t::Vector{Float64}, X::Matrix{Float64}, v::Vector{Float64})
  const tx = [t X]
  const (N,J) = size(tx)
  const m = sum(Int,v)
  const obs_tx = Matrix{Float64}(m,J)

  k = 0
  for i in 1:N
    if v[i] == 1
      k += 1
      for j in 1:J
        obs_tx[k,j] = tx[i,j]
      end
    end
  end

  return (obs_tx[:,1], obs_tx[:,2:end])
end

include("coxph_weibull.jl")
include("summary_cox.jl")

#=
#const Model = [:weibull, :piecewise, :gammaProcess]

function loglike_piecewise(t::Vector{Float64}, X::Matrix{Float64}, v::Vector{Float64})
end

function loglike_gammaProcess(t::Vector{Float64}, X::Matrix{Float64}, v::Vector{Float64})
end
=#

end # Cox
