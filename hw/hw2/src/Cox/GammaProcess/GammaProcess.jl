module GammaProcess

using Distributions, RCall
import MCMC, Base.show
# see slides 8 for more details on theory

export summary, plot, gp, Priorᵦ, Priorₕ, plotsurv, est_survival

include("gammaprocess.jl")
include("summary.jl")
include("plot.jl")
include("est_survival.jl")

end # GammaProcess
