println("Starting...")
include("Frailty/Frailty.jl")
using RCall
sym(M::Matrix{Float64}) = (M' + M) / 2
srand(276);

#=
model:
# h(tᵢⱼ| wᵢ,xᵢⱼ) = h₀(tᵢⱼ) wᵢ exp(xᵢⱼ'β)
=#
R"""
library(survival)
kidney <- read.table("dat/kidney.txt",header=TRUE,skip=5)
kidney$sex <- ifelse(kidney$sex==2, 0, 1) # M=1, F=0
cox_mod <- coxph(Surv(time,nu) ~ sex+age, data=kidney)
frail_mod <- coxph(Surv(time,nu) ~ sex+age+frailty(cluster,theta=1), data=kidney)
""";
@rget kidney;

const t = kidney[:time] * 1.;
const X = [kidney[:sex] kidney[:age]] * 1.; # sex, age
const v = kidney[:nu] * 1.;
const group = convert(Vector{Int64},kidney[:cluster]);

prior_β = Frailty.Prior_β(zeros(2), eye(2)*10^3, sym(inv(X'X))*10.)
prior_λ = Frailty.Prior_λ(.001,.001)
prior_α = Frailty.Prior_α(.001,.001,.1)
prior_η = Frailty.Prior_η(.001,.001,.1)

println(R"frail_mod")
@time out = Frailty.fit(t,X,v,group,prior_β,prior_λ,prior_α,prior_η,2000,10000)

β = hcat(map(m -> m.β, out)...)'
mean(β,1)
std(β,1)
size(unique(β,1),1) / size(β,1)

#=
include("proj1.jl")
=#
