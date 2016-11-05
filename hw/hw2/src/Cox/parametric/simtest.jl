include("Cox.jl")
using Distributions
using RCall
R"library(survival)";
srand(1);

function rweib(α::Float64, λ::Float64)
  assert(α > 0 && λ > 0)
  return (-log(1-rand())/λ)^(1/α)
end

function rweib(α::Float64, λ::Float64, n::Int)
  const out = Vector{Float64}(n)
  for i in 1:n
    out[i] = rweib(α,λ)
  end
  return out
end

const N = 1000
const β = [.2,.3,.4,.5]
const J = length(β)
const v = rand(Bernoulli(.6),N) * 1.
const X = randn(N,J)
const λ = 1.
const α = 2.
const t_uncensored = [rweib(α, exp(xb)*λ) for xb in X*β]
const t = [ v[i] == 1 ? t_uncensored[i] : rand() * t_uncensored[i] for i in 1:N]

@time m = Cox.coxph_weibull(t,X,v,
                            Cox.Diag(fill(1.,J+2)*1E-3),
                            B=20000,burn=10000);
s = Cox.summary(m)
Cox.plot(m);

@rput t X v;
println(R"summary(coxph(Surv(t,v) ~ X))")
