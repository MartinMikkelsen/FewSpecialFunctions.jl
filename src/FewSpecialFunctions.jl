module FewSpecialFunctions

export complex_quadrature
include("Complex_quad.jl")

export regular_coulomb, irregular_coulomb, C
include("CoulombWave.jl")

export Debye_function
include("Debye.jl")

export Fresnel_S_integral_pi, Fresnel_C_integral_pi , Fresnel_S_integral, Fresnel_C_integral , Fresnel_S_err, Fresnel_C_err
include("Fresnel.jl")

export Struve
include("Struve.jl")

export Clausen
include("Clausen.jl")
end
