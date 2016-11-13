module Frailty

using Distributions, RCall
import MCMC, Base.show, Base.summary

export plot, summary, est_survival, fit, Prior_β, Prior_λ, Prior_α, Prior_η

include("fit.jl")
include("summary.jl")
include("plot.jl")
include("est_survival.jl")

end #Frailty
