module FewSpecialFunctions

export complex_quadrature
include("Complex_quad.jl")

export regular_coulomb, irregular_coulomb, θ, Coulomb_H_minus, Coulomb_H_plus, Coulomb_cross, regular_Coulomb_approx, irregular_Coulomb_approx, regular_Coulomb_limit, irregular_Coulomb_limit
include("CoulombWave.jl")

export Debye_function
include("Debye.jl")

export Fresnel_S_integral_pi, Fresnel_C_integral_pi , Fresnel_S_integral, Fresnel_C_integral , Fresnel_S_erf, Fresnel_C_erf
include("Fresnel.jl")

export Struve
include("Struve.jl")

export Clausen
include("Clausen.jl")
end
