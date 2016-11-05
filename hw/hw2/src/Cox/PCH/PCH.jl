module PCH # piecewise-constant hazards

using Distributions, RCall
import MCMC.gibbs, Base.show, Base.summary

export summary, plot, pch, Priorᵦ, Priorᵧ

include("pch.jl")
include("summary.jl")
include("plot.jl")

end # PCH
