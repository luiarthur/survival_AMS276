immutable Summary_lm
  β̄::Vector{Float64}
  std_β::Vector{Float64}
  quantile_β::Matrix{Float64}

  ᾱ::Float64
  std_α::Float64
  quantile_α::Vector{Float64}

  λ̄::Float64
  std_λ::Float64
  quantile_λ::Vector{Float64}

  acc::Float64
  DIC::Float64
  B::Int
end

const quants = [0., .025, .25, .5, .75, .975, 1.]

function summary_cox(model::Cox_weib)
  const β= hcat(map(o -> o.β, model.params)...)'
  const (B,P) = size(β)
  const β̄= vec(mean(β,1))
  const std_β = vec(std(β,1))
  const quantile_β = hcat([ quantile(β[:,p], quants) for p in 1:P]...)

  const α = map(o -> sqrt(o.α), model.params)
  const ᾱ = mean(α)
  const std_α = std(α)
  const quantile_α = quantile(α, quants)

  const λ = map(o -> sqrt(o.λ), model.params)
  const λ̄ = mean(λ)
  const std_λ = std(λ)
  const quantile_λ = quantile(λ, quants)

  const acc = length(unique(α)) / B
  const DIC = 0.

  return Summary_lm(β̄,std_β,quantile_β,
                    ᾱ,std_α,quantile_α,
                    λ̄,std_λ,quantile_λ,
                    acc,DIC,B)
end

function show(io::IO, SL::Summary_lm)
  const P = length(SL.β̄)
  zeroInCI = [SL.quantile_β[2,p] <= 0 <= SL.quantile_β[6,p] ? "" : "*" for p in 1:P]

  @printf "%5s%10s%10s%3s\n" "" "mean" "std" "≠0"
  for k in 1:P
    @printf "%5s%10.4f%10.4f%3s\n" string("β",k) SL.β̄[k] SL.std_β[k] zeroInCI[k]
  end
  @printf "%5s%10.4f%10.4f\n" "α" SL.ᾱ SL.std_α
  @printf "%5s%10.4f%10.4f\n" "λ" SL.λ̄ SL.std_λ
  println(repeat("-",28))
  @printf "%5s%10.2f\n" "acc" SL.acc
  @printf "%5s%10.2f\n" "DIC" SL.DIC
end
