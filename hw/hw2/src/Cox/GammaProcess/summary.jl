const quants = [0., .025, .25, .5, .75, .975, 1.]

immutable Summary_GP
  β̄::Vector{Float64}
  std_β::Vector{Float64}
  quantile_β::Matrix{Float64}

  h̄::Vector{Float64}
  std_h::Vector{Float64}
  quantile_h::Matrix{Float64}

  accβ::Float64
  acch::Float64

  DIC::Float64

  B::Int
end

function summary(states::Vector{State})
  const β = hcat(map(s -> s.β, states)...)'
  const (B,P) = size(β)
  const β̄ = vec(mean(β,1))
  const std_β = vec(std(β,1))
  const quantile_β = hcat([ quantile(β[:,p], quants) for p in 1:P]...)'

  const h = hcat(map(s -> s.h, states)...)'
  const J = size(h, 2)
  const h̄ = vec(mean(h,1))
  const std_h = vec(std(h,1))
  const quantile_h = hcat([ quantile(h[:,j], quants) for j in 1:J]...)'

  const accβ = size(unique(β,1),1) / B
  const acch = size(unique(h,1),1) / B
  const DIC = 0.

  return Summary_GP(β̄, std_β, quantile_β,
                    h̄, std_h, quantile_h,
                    accβ, acch, DIC, B)
end

function show(io::IO, S::Summary_GP)
  const P = length(S.β̄)
  const J = length(S.h̄)
  zeroInCI = S.quantile_β[:,2] .<= 0 .<= S.quantile_β[:,6]

  @printf "%5s%10s%10s%10s%10s%3s\n" "" "mean" "std" "lower" "upper" "≠0"
  for p in 1:P
    @printf "%5s%10.4f%10.4f%10.4f%10.4f%3s\n" string("β",p) S.β̄[p] S.std_β[p] S.quantile_β[p,2] S.quantile_β[p,6] zeroInCI[p] ? "" : "*"
  end

  println(repeat("-",48))
  for j in 1:J
    @printf "%5s%10.4f%10.4f%10.4f%10.4f\n" string("h",j) S.h̄[j] S.std_h[j] S.quantile_h[j,2] S.quantile_h[j,6]
  end

  println(repeat("-",48))
  @printf "%5s%10.2f\n" "accβ:" S.accβ
  @printf "%5s%10.2f\n" "acch:" S.acch
  @printf "%5s%10.2f\n" "DIC" S.DIC
end
