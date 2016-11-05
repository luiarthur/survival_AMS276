module SurvivalDistributions

import Distributions: Normal, Logistic, logpdf, logccdf

export logpdf, logccdf, lognormal, loglogistic, weibull

immutable ExtremeValue ; end
logpdf{T <: Real}(E::ExtremeValue, z::T) = z - exp(z)
logccdf{T <: Real}(E::ExtremeValue, z::T) = -exp(z)
ErrorDistribution = Union{Normal, Logistic, ExtremeValue}


immutable TimeDistribution
  errorDist::ErrorDistribution
end
lognormal = TimeDistribution(Normal())
loglogistic = TimeDistribution(Logistic())
weibull = TimeDistribution(ExtremeValue())


function logpdf{T <: TimeDistribution}(td::T, t::Float64, m::Float64, s::Float64)
  assert(t > 0 && s > 0)
  logpdf(td.errorDist, (log(t)-m)/s) - log(t*s)
end

function logccdf{T <: TimeDistribution}(td::T, t::Float64, m::Float64, s::Float64)
  assert(t > 0 && s > 0)
  logccdf(td.errorDist, (log(t)-m) / s)
end

end # module SurvivalDistributions

#= 
include("SurvivalDistributions.jl")
=# 
