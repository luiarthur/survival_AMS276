module PCH # piecewise-constant hazards

using Distributions, RCall
import MCMC, Base.show

export summary, plot, pch, Priorβ, Priorλ, plotsurv

include("pch.jl")
include("summary.jl")
include("plot.jl")
include("plotsurv.jl")

end # PCH
