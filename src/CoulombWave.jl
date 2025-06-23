using SpecialFunctions
using HypergeometricFunctions

"""
    η(a::Number, k::Number)
    η(ϵ::Number)

Coulomb parameter. For two arguments, returns 1/(a*k). For one argument, returns 1/sqrt(ϵ).
"""
function η(a::Number, k::Number)
    return 1/(a*k)
end

function η(ϵ::Number)
    return 1/sqrt(ϵ)
end

"""
    C(ℓ::Number, η::Number)

Coulomb normalization constant.
"""
function C(ℓ::Number, η::Number)
    logg = loggamma(ℓ+1+im*η)
    return 2^ℓ * exp(-π*η/2) * exp(real(logg)) / gamma(2*ℓ+2)
end

"""
    θ(ℓ::Number, η::Number, ρ::Number)

Coulomb phase function.
"""
function θ(ℓ::Number, η::Number, ρ::Number)
    logg = loggamma(ℓ+1+im*η)
    return ρ - ℓ*π*0.5 - η*log(2*ρ) + imag(logg)
end

"""
    F(ℓ::Number, η::Number, ρ::Number)

Regular Coulomb wave function.

References:
- [Coulomb wave function](https://en.wikipedia.org/wiki/Coulomb_wave_function)
- [Implementation paper](https://arxiv.org/abs/1804.10976)
"""
function F(ℓ::Number, η::Number, ρ::Number)
    ℓc = complex(float(ℓ))
    ηc = complex(float(η))
    ρc = complex(float(ρ))
    return C(ℓc, ηc) * ρc^(ℓc+1) * exp(-im*ρc) * _₁F₁(complex(ℓc+1-im*ηc), complex(2*ℓc+2), complex(2*im*ρc))
end

function F(ℓ::Real, η::Real, ρ::Real)
    return real(F(complex(ℓ), complex(η), complex(ρ)))
end

"""
    D⁺(ℓ::Number, η::Number)

Coulomb D⁺ normalization factor.
"""
function D⁺(ℓ::Number, η::Number)
    return (-2*im)^(2*ℓ+1) * gamma(ℓ+1+im*η) / (C(ℓ,η) * gamma(2*ℓ+2))
end

"""
    D⁻(ℓ::Number, η::Number)

Coulomb D⁻ normalization factor.
"""
function D⁻(ℓ::Number, η::Number)
    return (2*im)^(2*ℓ+1) * gamma(ℓ+1-im*η) / (C(ℓ,η) * gamma(2*ℓ+2))
end

"""
    H⁺(ℓ::Number, η::Number, ρ::Number)

Outgoing Coulomb wave function.

References:
- [Coulomb wave function](https://en.wikipedia.org/wiki/Coulomb_wave_function)
- [Implementation paper](https://arxiv.org/abs/1804.10976)
"""
function H⁺(ℓ::Number, η::Number, ρ::Number)
    return D⁺(ℓ,η) * ρ^(ℓ+1) * exp(+im*ρ) * HypergeometricFunctions.U(ℓ+1+im*η, 2*ℓ+2, -2*im*ρ)
end

"""
    H⁻(ℓ::Number, η::Number, ρ::Number)

Incoming Coulomb wave function.

References:
- [Coulomb wave function](https://en.wikipedia.org/wiki/Coulomb_wave_function)
- [Implementation paper](https://arxiv.org/abs/1804.10976)
"""
function H⁻(ℓ::Number, η::Number, ρ::Number)
    return D⁻(ℓ,η) * ρ^(ℓ+1) * exp(-im*ρ) * HypergeometricFunctions.U(ℓ+1-im*η, 2*ℓ+2, +2*im*ρ)
end

"""
    F_imag(ℓ::Number, η::Number, ρ::Number)

Imaginary part of the regular Coulomb wave function.
"""
function F_imag(ℓ::Number, η::Number, ρ::Number)
    return (H⁺(ℓ, η, ρ) - H⁻(ℓ, η, ρ)) / (2*im)
end

"""
    G(ℓ::Number, η::Number, ρ::Number)

Irregular Coulomb wave function.

References:
- [Coulomb wave function](https://en.wikipedia.org/wiki/Coulomb_wave_function)
- [Implementation paper](https://arxiv.org/abs/1804.10976)
"""
function G(ℓ::Number, η::Number, ρ::Number)
    return (H⁺(ℓ, η, ρ) + H⁻(ℓ, η, ρ)) / 2
end

# Explicit real-valued overload to guarantee real output and avoid complex issues in downstream code
function G(ℓ::Real, η::Real, ρ::Real)
    return real(G(complex(ℓ), complex(η), complex(ρ)))
end

"""
    M_regularized(α::Number, β::Number, γ::Number)

Regularized confluent hypergeometric function.
"""
function M_regularized(α::Number, β::Number, γ::Number)
    return 1/gamma(β) * HypergeometricFunctions._₁F₁(α, β, γ)
end

"""
    Φ(ℓ::Number, η::Number, ρ::Number)

Modified Coulomb function Φ.
"""
function Φ(ℓ::Number, η::Number, ρ::Number)
    return (2*η*ρ)^(ℓ+1) * exp(im*ρ) * M_regularized(ℓ+1+im*η, 2*ℓ+2, -2*im*ρ)
end

"""
    w(ℓ::Integer, η::Number)
    w(ℓ::Number, η::Number)

Auxiliary function for Coulomb wave functions.
"""
function w(ℓ::Integer, η::Number)
    result = one(η)
    for j in 0:ℓ
        result *= 1 + j^2/η^2
    end
    return result
end

function w(ℓ::Number, η::Number)
    if isapprox(mod(ℓ-0.5, 1.0), 0.0; atol=1e-12)
        result = one(η)
        j = 0.5
        while j <= ℓ
            result *= 1 + j^2/η^2
            j += 1
        end
        return result
    else
        throw(ArgumentError("ℓ must be either an integer or half-integer (1/2, 3/2, 5/2, ...)"))
    end
end

"""
    w_plus(ℓ::Number, η::Number)

Auxiliary function for Coulomb wave functions.
"""
function w_plus(ℓ::Number, η::Number)
    return gamma(ℓ+1+im*η)/( (im*η)^(2*ℓ+1) * gamma(-ℓ+im*η) )
end

"""
    w_minus(ℓ::Number, η::Number)

Auxiliary function for Coulomb wave functions.
"""
function w_minus(ℓ::Number, η::Number)
    return gamma(ℓ+1-im*η)/( (-im*η)^(2*ℓ+1) * gamma(-ℓ-im*η) )
end

"""
    h_plus(ℓ::Number, η::Number)

Auxiliary function for Coulomb wave functions.
"""
function h_plus(ℓ::Number, η::Number)
    return (digamma(ℓ+1+im*η) + digamma(-ℓ+im*η))/2 - log(im*η)
end

"""
    h_minus(ℓ::Number, η::Number)

Auxiliary function for Coulomb wave functions.
"""
function h_minus(ℓ::Number, η::Number)
    return (digamma(ℓ+1-im*η) + digamma(-ℓ-im*η))/2 - log(-im*η)
end

"""
    g(ℓ::Number, η::Number)

Auxiliary function for Coulomb wave functions.
"""
function g(ℓ::Number, η::Number)
    x = (digamma(ℓ+1+im*η) + digamma(ℓ+1-im*η)) / 2 - log(abs(η))
    return real(x)
end

"""
    Φ_dot(ℓ::Number, η::Number, ρ::Number; h=1e-6)

Numerical derivative of Φ with respect to ℓ.
"""
function Φ_dot(ℓ::Number, η::Number, ρ::Number; h=1e-6)
    return (Φ(ℓ+h, η, ρ) - Φ(ℓ-h, η, ρ))/(2h)
end

"""
    F_dot(ℓ::Number, η::Number, ρ::Number; h=1e-6)

Numerical derivative of F with respect to ℓ.
"""
function F_dot(ℓ::Number, η::Number, ρ::Number; h=1e-6)
    return (F(ℓ+h, η, ρ) - F(ℓ-h, η, ρ))/(2h)
end

"""
    Ψ(ℓ::Number, η::Number, ρ::Number; h=1e-6)

Auxiliary function for Coulomb wave functions.
"""
function Ψ(ℓ::Number, η::Number, ρ::Number; h=1e-6)
    return w(ℓ,η)*Φ_dot(ℓ,η,ρ;h=h)/2 + Φ_dot(-ℓ-1,η,ρ;h=h)/2
end

"""
    I(ℓ::Number, η::Number, ρ::Number; h=1e-6)

Auxiliary function for Coulomb wave functions.
"""
function I(ℓ::Number, η::Number, ρ::Number; h=1e-6)
    return C(ℓ,η)*gamma(2*ℓ+2)/((2*η)^(ℓ+1)) * Ψ(ℓ,η,ρ;h=h)
end


export η, C, θ, F, D⁺, D⁻, H⁺, H⁻, F_imag, G, M_regularized, Φ, w, w_plus, w_minus, h_plus, h_minus, g, Φ_dot, F_dot, Ψ, I