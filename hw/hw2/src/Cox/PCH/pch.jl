immutable Priorᵦ
  mean::Vector{Float64}
  cov::Matrix{Float64}
  proposal::Matrix{Float64}
end

immutable Priorᵧ
  mean::Vector{Float64}
  var::Vector{Float64}
end

#function pch(y::Vector{Float64}, X::Matrix{Float64}, v::Vector{Float64}, 
#             grid::Vector{Float64},)


