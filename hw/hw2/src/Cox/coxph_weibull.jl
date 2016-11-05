immutable Param_weib
  α::Float64
  λ::Float64
  β::Vector{Float64}
end

immutable Cox_weib
  params::Vector{Param_weib}
  t::Vector{Float64}
  X::Matrix{Float64}
  v::Vector{Float64}
end

function loglike_weibull(t::Vector{Float64}, X::Matrix{Float64}, v::Vector{Float64})
  #const (t_obs, X_obs) = part(t,X,v)
  #log_h(p::Param_weib) = sum(p.α * log(t_obs) + log(p.α*p.λ) + X_obs*p.β)
  log_h(p::Param_weib) = sum(v .* (p.α * log(t) + log(p.α*p.λ) + X*p.β))
  log_S(p::Param_weib) = -sum(exp(X * p.β) .* t.^p.α .* p.λ)
  ll(p::Param_weib) = log_h(p) + log_S(p)

  return ll
end

logprior_β(β::Vector{Float64}) = 0. #(- β'β / 200)[1]
logprior_α(α::Float64) = logpdf(Gamma(.1,10),α) #logpdf(Gamma(.1,10),α)
logprior_λ(λ::Float64) = logpdf(Gamma(.1,10),λ)

function coxph_weibull(t::Vector{Float64}, X::Matrix{Float64}, v::Vector{Float64},
                       Σ::Matrix{Float64}; B::Int=2000, burn::Int=100)

  const J = size(X,2)

  function logprior_weibull(p::Param_weib)
    return logprior_β(p.β) + logprior_α(p.α) + logprior_λ(p.λ)
  end
  ll = loglike_weibull(t,X,v)

  ll_plus_lp(p::Param_weib) = logprior_weibull(p) + ll(p)

  const lower_bound = [ 0.;  0.; fill(-Inf,J)]
  const upper_bound = [Inf; Inf; fill( Inf,J)]

  function update(curr_state::Param_weib)
    const curr = [curr_state.α; curr_state.λ; copy(curr_state.β)]
    const cand = rand(MvNormal(curr, Σ))
    const cand_state = Param_weib(cand[1], cand[2], cand[3:end])

    const in_bounds = all(lower_bound .< cand .< upper_bound)

    if in_bounds && log(rand()) < ll_plus_lp(cand_state) - ll_plus_lp(curr_state)
      const new_state = cand_state
    else
      const new_state = curr_state
    end

    return new_state
  end

  const init = Param_weib(1., 1., zeros(Float64,J))
  return Cox_weib(gibbs(init,update,B,burn),t,X,v)
end # cox_weibull
