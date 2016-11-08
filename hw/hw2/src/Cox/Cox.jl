module Cox

export Diag

Diag{T <: Real}(v::Vector{T}) = Matrix(Diagonal(v))

include("Parametric/Parametric.jl"); using .Parametric
include("PCH/PCH.jl"); using .PCH
include("GammaProcess/GammaProcess.jl"); using .GammaProcess

end # Cox
