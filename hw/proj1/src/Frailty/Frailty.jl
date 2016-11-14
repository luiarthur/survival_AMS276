module Frailty

using Distributions, RCall
import MCMC, Base.show

export plot, summary, est_survival, fit, Prior_β, Prior_λ, Prior_α, Prior_η

R"""
# To install rcommon:
# install.packages("devtools")
# devtools::install_github("luiarthur/rcommon")
library("rcommon")
"""

include("fit.jl")
include("summary.jl")
include("plot.jl")
include("est_survival.jl")

end #Frailty
