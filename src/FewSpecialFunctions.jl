module FewSpecialFunctions

export η, C, θ, F, D⁺, D⁻, H⁺, H⁻, F_imag, G, M_regularized, Φ, w
include("CoulombWave.jl")

export debye_function
include("Debye.jl")

export Fresnel_S_integral_pi, Fresnel_C_integral_pi , Fresnel_S_integral, Fresnel_C_integral , Fresnel_S_erf, Fresnel_C_erf
include("Fresnel.jl")

export Struve
include("Struve.jl")

export Clausen
include("Clausen.jl")

export FermiDiracIntegral, FermiDiracIntegralNorm
include("FermiDirac.jl")

export MarcumQ, dQdb
include("MarcumQ.jl")

export U, V, W, dU, dV, dW
include("parabolic_cylinder.jl")
end
