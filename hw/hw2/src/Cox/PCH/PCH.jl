module PCH # piecewise-constant hazards

using Distributions, RCall
import MCMC, Base.show

# see slides7.pdf for the theory.

export summary, plot, pch, Priorβ, Priorλ, plotsurv, est_survival

include("pch.jl")
include("summary.jl")
include("plot.jl")
include("plotsurv.jl")
include("est_survival.jl")
include("plotCI.jl")

end # PCH
