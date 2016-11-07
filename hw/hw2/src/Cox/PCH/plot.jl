function plot(states::Vector{State})
  R"""
  # To install rcommon:
  # install.packages("devtools")
  # devtools::install_github("luiarthur/rcommon")
  library(rcommon)
  """

  beta = hcat(map(p -> p.β, states)...)'
  lambda = hcat(map(p -> p.λ, states)...)'

  @rput beta lambda

  R"""plotPosts(cbind(beta,lambda),show.x=FALSE,
  cnames=c(paste0('b',1:ncol(beta)),paste0('l',1:ncol(lambda))),stats=FALSE)
  """;
end
