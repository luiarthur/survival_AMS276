const quants = [0., .025, .25, .5, .75, .975, 1.]

immutable Summary_PCH
  β̄::Vector{Float64}
  std_β::Vector{Float64}
  quantile_β::Matrix{Float64}

  λ̄::Vector{Float64}
  std_λ::Vector{Float64}
  quantile_λ::Matrix{Float64}

  accβ::Float64
  accλ::Float64

  DIC::Float64

  B::Int
end

function summary(states::Vector{State})
  const β = hcat(map(s -> s.β, states)...)'
  const (B,P) = size(β)
  const β̄ = vec(mean(β,1))
  const std_β = vec(std(β,1))
  const quantile_β = hcat([ quantile(β[:,p], quants) for p in 1:P]...)'

  const λ = hcat(map(s -> s.λ, states)...)'
  const J = size(λ, 2)
  const λ̄ = vec(mean(λ,1))
  const std_λ = vec(std(λ,1))
  const quantile_λ = hcat([ quantile(λ[:,j], quants) for j in 1:J]...)'

  const accβ = size(unique(β,1),1) / B
  const accλ = size(unique(λ,1),1) / B
  const DIC = 0.

  return Summary_PCH(β̄, std_β, quantile_β,
                     λ̄, std_λ, quantile_λ,
                     accβ, accλ, DIC, B)
end

function show(io::IO, S::Summary_PCH)
  const P = length(S.β̄)
  const J = length(S.λ̄)
  zeroInCI = S.quantile_β[:,2] .<= 0 .<= S.quantile_β[:,6]

  @printf "%5s%10s%10s%10s%10s%3s\n" "" "mean" "std" "lower" "upper" "≠0"
  for p in 1:P
    @printf "%5s%10.4f%10.4f%10.4f%10.4f%3s\n" string("β",p) S.β̄[p] S.std_β[p] S.quantile_β[p,2] S.quantile_λ[p,6] zeroInCI[p] ? "" : "*"
  end

  println(repeat("-",48))
  for j in 1:J
    #@printf "%5s%10.4f%10.4f\n" string("λ",j) S.λ̄[j] S.std_λ[j]
    @printf "%5s%10.4f%10.4f%10.4f%10.4f\n" string("λ",j) S.λ̄[j] S.std_λ[j] S.quantile_λ[j,2] S.quantile_λ[j,6]
  end

  println(repeat("-",48))
  @printf "%5s%10.2f\n" "accβ:" S.accβ
  @printf "%5s%10.2f\n" "accλ:" S.accλ
  @printf "%5s%10.2f\n" "DIC" S.DIC
end
