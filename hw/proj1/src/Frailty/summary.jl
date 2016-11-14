function summary_v(v::Vector{Float64}; a::Float64=.05)
  const MEAN = mean(v)
  const SD = std(v)
  const ACC = length(unique(v)) / length(v)
  const q = quantile(v, [a/2, 1-a/2])

  return Summary_scalar(MEAN,SD,ACC,q)
end

function summary_vv(vv::Vector{Vector{Float64}}; a::Float64=.05)
  const J = length(vv[1])
  const m = hcat(vv...)'
  const MEAN = vec(mean(m,1))
  const SD = vec(std(m,1))
  const ACC = size(unique(m,1),1) / length(vv)
  const q = hcat([ quantile(m[:,j], [a/2, 1-a/2]) for j in 1:J ]...)'

  return Summary_vec(MEAN,SD,ACC,q)
end

immutable Summary_vec
  MEAN::Vector{Float64}
  SD::Vector{Float64}
  ACC::Float64
  q::Matrix{Float64}
end

immutable Summary_scalar
  MEAN::Float64
  SD::Float64
  ACC::Float64
  q::Vector{Float64}
end

immutable Summary_frailty
  β::Summary_vec
  λ::Summary_scalar
  α::Summary_scalar
  η::Summary_scalar
  w::Summary_vec
end

function summary(S::Vector{State}; a::Float64=.05)
  const summary_β = summary_vv(map(s->s.β, S), a=a)
  const summary_λ = summary_v(map(s->s.λ, S), a=a)
  const summary_α = summary_v(map(s->s.α, S), a=a)
  const summary_w = summary_vv(map(s->s.w, S), a=a)
  const summary_η = summary_v(map(s->s.η, S), a=a)
  
  return Summary_frailty(summary_β,summary_λ,summary_α,summary_η,summary_w)
end
