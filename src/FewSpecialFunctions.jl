module FewSpecialFunctions

export complex_quadrature
include("Complex_quad.jl")

export regular_coulomb, irregular_coulomb, C
include("CoulombWave.jl")

export Debye_function
include("Debye.jl")

export S_integral_pi, C_integral_pi ,S_integral, C_integral, S_err, C_err
include("Fresnel.jl")

export Struve
include("Struve.jl")
end
