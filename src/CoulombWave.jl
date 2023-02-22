export complex_quadrature, regular_coulomb, irregular_coulomb, C, θ, Coulomb_H_minus, Coulomb_H_plus, Coulomb_cross, regular_Coulomb_approx, irregular_Coulomb_approx, regular_Coulomb_limit, irregular_Coulomb_limit

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

@doc raw"""
    C(ℓ,η)

Returns Coulomb normalization constant given by
```math
    C_\ell(\eta) = \frac{2^\ell \exp{-\pi \eta/2} |\Gamma(\ell+1+i \eta)|}{(2\ell+1)!}
```
"""
function C(ℓ,η)
    return 2^ℓ*exp(-π*η/2).*(abs(gamma(ℓ+1+1im*η))/(factorial(2*ℓ+1)))
end
@doc raw"""
    θ(ℓ,η,ρ)

Returns the phase of the Coulomb functions given by
```math
    \theta_\ell(\eta,\rho) = \rho - \eta \ln(2\rho) - \frac{1}{2}\ell \pi + \sigma_\ell(\eta)
```
"""
function θ(ℓ,η,ρ)
    return ρ - η.*log(2*ρ)-0.5*ℓ*π+angle.(gamma(ℓ+1+1im*η))
end
@doc raw"""
    Coulomb_H_minus(ℓ,η,ρ)

Complex Coulomb wave function. Infinity handled using the substitution f(t) -> f(u/(1-u)*1/(1-u)^2).
Returns Coulomb wave function 
```math
    H^{-}_\ell = G_\ell - iF_\ell
```
"""
function Coulomb_H_minus(ℓ,η,ρ)
    return (exp(-1im.*ρ)*ρ.^(-ℓ))/(factorial(2*ℓ+1).*C(ℓ,η)).*complex_quadrature(t -> exp(-t)*t.^(ℓ-1im*η)*(t+2*1im*ρ).^(ℓ+1im*η),0,Inf)
end 
@doc raw"""
    irregular_coulomb(ℓ,η,ρ)

Irregular Coulomb wave function ℓ is the order(non-negative integer), η is the charge (real parameter) and ρ is the radial coordinate (non-negative real variable).

returns the value G_ℓ(η,ρ)
"""
function irregular_Coulomb(ℓ,η,ρ)
    return real(Coulomb_H_minus(ℓ,η,ρ))
end
@doc raw"""

Returns Coulomb wave function 
```math
    H^{+}_\ell = G_\ell + iF_\ell
```
"""
function Coulomb_H_plus(ℓ,η,ρ)
    return irregular_Coulomb(ℓ,η,ρ) .+ 1im*regular_coulomb(ℓ,η,ρ)
end
@doc raw"""
    Coulomb_cross(ℓ,η)

Wronskian relation / cross product.

```math
    F_{\ell-1}G_{\ell}-F_{\ell}G_{\ell-1} = \ell/(\ell^2+\eta^2)^{1/2}
```
"""
function Coulomb_cross(ℓ,η)
    if ℓ >= 1
        return ℓ/(ℓ^2+η^2)^(0.5)
    end
    if ℓ < 1
        return println("ℓ must be larger than or equal to 1")
    end
end
@doc raw"""
    regular_Coulomb_approx(ℓ,η,ρ)

For ρ -> 0 and η fixed approximate the regular Coulomb wave function as
```math
    F_\ell(\eta,\rho) \simeq C_\ell(\eta)^{\ell+1}
```
"""
function regular_Coulomb_approx(ℓ,η,ρ)
    return C(ℓ,η)*ρ^(ℓ+1)
end 
@doc raw"""
    irregular_Coulomb_approx(ℓ,η,ρ)

For ρ -> 0 and η fixed approximate the irregular Coulomb wave function as
```math
    G_\ell(\eta,\rho) \simeq \frac{\rho^{-\ell}}{(2\ell+1)C_\ell(\eta)}
```
"""
function irregular_Coulomb_approx(ℓ,η,ρ)
    return ρ^(-ℓ)/((2*ℓ+1)*C(ℓ,η))
end
@doc raw"""
    regular_Coulomb_limit(ℓ,η,ρ)

In the limit ρ -> ∞ with η fixed, returns the regular Coulomb wave as 
```math
    F_{\ell}(\eta,\rho) \simeq \sin(\theta_\ell(\eta,\rho))
```
"""
function regular_Coulomb_limit(ℓ,η,ρ)
    return sin.(θ(ℓ,η,ρ))
end 
@doc raw"""
    irregular_Coulomb_limit(ℓ,η,ρ)

In the limit ρ -> ∞ with η fixed, returns the irregular Coulomb wave as
```math
    G_{\ell}(\eta,\rho) \simeq \cos(\theta_\ell(\eta,\rho))
```
"""
function irregular_Coulomb_limit(ℓ,η,ρ)
    return cos.(θ(ℓ,η,ρ))
end 