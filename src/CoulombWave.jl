export complex_quadrature, regular_coulomb, irregular_coulomb, C

using SpecialFunctions #For gamma function

@doc raw"""
    regular_coulomb(ℓ,η,ρ)

Regular Coulomb wave function ℓ is the order(non-negative integer), η is the charge (real parameter) and ρ is the radial coordinate (non-negative real variable).

returns the value F_ℓ(η,ρ)
"""
function regular_coulomb(ℓ,η,ρ)
    First = ρ.^(ℓ+1)*2^ℓ*exp(1im.*ρ-(π.*η/2))/(abs(gamma(ℓ+1+1im*η)))
    Integral_value = complex_quadrature(t -> exp(-2*1im.*ρ*t).*t.^(ℓ+1im*η)*(1-t).^(ℓ-1im*η),1e-5,1)
    return real(First.*Integral_value)
end

#Needs work
function C(ℓ::Int,η::Float64)
    return 2^ℓ*exp(-π*η/2)*(abs(gamma(ℓ+1+1im*η))/(factorial(2*ℓ+1)))
end

function irregular_coulomb(ℓ,η,ρ)
    First = exp(-1im*ρ)*ρ^(-ℓ)/(factorial(2*ℓ+1).*C(ℓ,η))
    Integral_value = complex_quadrature(t -> exp(-t)*t^(-ℓ-1im*η)*(t+2*1im*ρ)^(ℓ+1im*η),1e-3,1e5)
    print(real(First*Integral_value))
    return First.*Integral_value
end
