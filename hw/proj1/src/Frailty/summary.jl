function summary_v(v::Vector{Float64}; a::Float64=.05)
  const MEAN = mean(v)
  const SD = std(v)
  const ACC = length(unique(v)) / length(v)
  const q = quantile(v, [a/2, 1-a/2])

  return (MEAN,SD,ACC,q)
end

function summary_vv(vv::Vector{Vector{Float64}}; a::Float64=.05)
  const J = length(vv[1])
  const m = hcat(vv...)'
  const MEAN = mean(m,1)
  const SD = std(m,1)
  const ACC = size(unique(m,1),1) / length(vv)
  const q = hcat([ quantile(m[:,j], [a/2, 1-a/2]) for j in 1:J ]...)'

  return (MEAN,SD,ACC,q)
end

#immutable Summary_frailty
#  β::
#end


function summary(S::Vector{State}; a::Float64=.05)
  const summary_β = summary(map(m->m.β, S), a)
  const summary_λ = summary(map(m->m.λ, S), a)
  const summary_α = summary(map(m->m.α, S), a)
  const summary_w = summary(map(m->m.w, S), a)
  const summary_η = summary(map(m->m.η, S), a)
  
  # return Summary_frailty()
end
