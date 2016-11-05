module PCH # piecewise-constant hazards

using Distributions, RCall
import MCMC, Base.show, Base.summary

export summary, plot, pch, Priorβ, Priorλ

include("pch.jl")
include("summary.jl")
include("plot.jl")

end # PCH
