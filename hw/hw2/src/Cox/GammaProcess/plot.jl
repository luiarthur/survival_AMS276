function load_rcommon()
  R"""
  # To install rcommon:
  # install.packages("devtools")
  # devtools::install_github("luiarthur/rcommon")
  library(rcommon)
  """
end

function plot(states::Vector{State})
  load_rcommon()

  beta = hcat(map(p -> p.β, states)...)'
  h = hcat(map(p -> p.h, states)...)'

  @rput beta h

  R"""plotPosts(cbind(beta,h),show.x=FALSE,
  cnames=c(paste0('b',1:ncol(beta)),paste0('l',1:ncol(h))),stats=FALSE)
  """;
end

function plot(states::Vector{State}, param::String, k::Int)
  assert(in(param,["beta","h"]))

  load_rcommon()

  const v = param == "beta" ?  map(p->p.β[k],states) : map(p->p.h[k],states)

  @rput v k param
  R"plotPost(v,main=paste0(param,k))"
end

function plot(states::Vector{State}, param::String, k::Vector{Int})
  assert(in(param,["beta","h"]))

  load_rcommon()

  if param == "beta"
    m = hcat(map(p->p.β[k],states)...)'
  else
    m = hcat(map(p->p.h[k],states)...)'
  end

  @rput m k param
  R"plotPosts(m,cnames=paste0(param,k))"
end
