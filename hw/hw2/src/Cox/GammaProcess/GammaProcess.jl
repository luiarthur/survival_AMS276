module GammaProcess

using Distributions, RCall
import MCMC, Base.show

export summary, plot, pch, Priorβ, Priorλ, plotsurv, est_survival

include("gammaprocess.jl")
include("summary.jl")
include("plot.jl")
include("plotsurv.jl")
include("est_survival.jl")

end # GammaProcess
