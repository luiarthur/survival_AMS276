module Parametric

using Distributions, RCall
import MCMC.gibbs, Base.show, Base.summary

export coxph_weibull, summary, plot

include("coxph_weibull.jl")
include("summary.jl")
include("plot.jl")

end
