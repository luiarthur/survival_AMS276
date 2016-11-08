function plot(cox::Cox_weib)
  R"""
  # To install rcommon:
  # install.packages("devtools")
  # devtools::install_github("luiarthur/rcommon")
  library(rcommon)
  """

  params = cox.params
  beta = hcat(map(p -> p.β, params)...)'
  alpha = map(p -> p.α, params)
  lambda = map(p -> p.λ, params)

  @rput beta alpha lambda

  R"plotPosts(cbind(beta,alpha,lambda),show.x=FALSE)";
end

function est_survival(cox::Cox_weib,t::Vector{Float64},
                      xs::Vector{Vector{Float64}}, f)
  const S = est_survival(cox,t,xs)
  return mapslices(f, S, 1)[1,:,:]
end

function est_survival(cox::Cox_weib,t::Vector{Float64},xs::Vector{Vector{Float64}})
  return cat(3,[est_survival(cox,t,x0) for x0 in xs]...)
end

function est_survival(cox::Cox_weib,t::Vector{Float64},x0::Vector{Float64})
  params = cox.params
  β = hcat(map(p -> p.β, params)...)'
  α = map(p -> p.α, params)
  λ = map(p -> p.λ, params)
  const B = length(params)
  const N = length(t)

  return [ exp(-λ[b] * t[i]^α[b] * exp((x0'β[b])[1]) ) for b in 1:B, i in 1:N]
end

plotsurv = begin
  R"""
  plotsurv_parametric <- function(grid,S,col_l=rep('black',ncol(S)),add=FALSE,...) {
    stopifnot(nrow(S) == length(grid))
    stopifnot(min(S) >= 0 && max(S) <= 1)
    if (!add) plot(grid,xlim=c(-.1,max(grid)),ylim=c(0,1),type='n',...)
    P <- NCOL(S)
    J <- NROW(S)
    for (p in 1:P) lines(grid,S[,p],col=col_l[p],...)
  }
  """
  R"plotsurv_parametric"
end
