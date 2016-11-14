# REUSEABLE ###########################################

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

function pretty_print(s::Summary_scalar,name::String)
  @printf "%5s%10.4f%10.4f%10.4f%10.4f\n" name s.MEAN s.SD s.q[1] s.q[2]
end

# SPECIFIC ###########################################
dash(n::Int=48) =   println(repeat("-",48))
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

function show(io::IO,s::Summary_frailty)
  const P = length(s.β.MEAN)
  const N = length(s.w.MEAN)
  const zeroInCI = s.β.q[:,1] .<= 0 .<= s.β.q[:,2]

  @printf "%5s%10s%10s%10s%10s%3s\n" "" "mean" "std" "lower" "upper" "≠0"
  #println("------   -------   -------   -------   -------")
  for p in 1:P
    @printf "%5s%10.4f%10.4f%10.4f%10.4f%3s\n" string("β",p) s.β.MEAN[p] s.β.SD[p] s.β.q[p,1] s.β.q[p,2] zeroInCI[p] ? "" : "*"
  end
  dash()

  pretty_print(s.λ,"λ")
  pretty_print(s.α,"α")
  pretty_print(s.η,"η")

  dash()
  @printf "%5s%10.4f\n" "accβ" s.β.ACC
  @printf "%5s%10.4f\n" "accα" s.α.ACC
  @printf "%5s%10.4f\n" "accη" s.η.ACC
  dash()

  for i in 1:N
    @printf "%5s%10.4f%10.4f%10.4f%10.4f\n" string("w",i) s.w.MEAN[i] s.w.SD[i] s.w.q[i,1] s.w.q[i,2]
  end


  dash()
end
