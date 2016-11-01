function plot(cox::Cox_weib)
  params = cox.params
  beta = hcat(map(p -> p.β, params)...)'
  alpha = map(p -> p.α, params)
  lambda = map(p -> p.λ, params)

  @rput beta alpha lambda

  R"plotPosts(cbind(beta,alpha,lambda),show.x=FALSE)";
end
