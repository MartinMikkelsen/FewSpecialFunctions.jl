using SpecialFunctions
using HypergeometricFunctions

function η(a::Number, k::Number)
    return 1/(a*k)
end

function η(ϵ::Number)
    return 1/sqrt(ϵ)
end

function C(ℓ::Number, η::Number)
    logg = loggamma(ℓ+1+im*η)
    return 2^ℓ * exp(-π*η/2) * exp(real(logg)) / gamma(2*ℓ+2)
end

function θ(ℓ::Number, η::Number, ρ::Number)
    logg = loggamma(ℓ+1+im*η)
    return ρ - ℓ*π*0.5 - η*log(2*ρ) + imag(logg)
end

function F(ℓ::Number, η::Number, ρ::Number)
    return C(ℓ, η) * ρ^(ℓ+1) * exp(-im*ρ) * _₁F₁(ℓ+1-im*η, 2*ℓ+2, 2*im*ρ)
end

function D⁺(ℓ::Number, η::Number)
    return (-2*im)^(2*ℓ+1) * gamma(ℓ+1+im*η) / (C(ℓ,η) * gamma(2*ℓ+2))
end

function D⁻(ℓ::Number, η::Number)
    return (2*im)^(2*ℓ+1) * gamma(ℓ+1-im*η) / (C(ℓ,η) * gamma(2*ℓ+2))
end

function H⁺(ℓ::Number, η::Number, ρ::Number)
    return D⁺(ℓ,η) * ρ^(ℓ+1) * exp(+im*ρ) * HypergeometricFunctions.U(ℓ+1+im*η, 2*ℓ+2, -2*im*ρ)
end

function H⁻(ℓ::Number, η::Number, ρ::Number)
    return D⁻(ℓ,η) * ρ^(ℓ+1) * exp(-im*ρ) * HypergeometricFunctions.U(ℓ+1-im*η, 2*ℓ+2, +2*im*ρ)
end

function F_imag(ℓ::Number, η::Number, ρ::Number)
    return (H⁺(ℓ, η, ρ) - H⁻(ℓ, η, ρ)) / (2*im)
end

function G(ℓ::Number, η::Number, ρ::Number)
    return (H⁺(ℓ, η, ρ) + H⁻(ℓ, η, ρ)) / 2
end

function M_regularized(α::Number, β::Number, γ::Number)
    return 1/gamma(β) * HypergeometricFunctions._₁F₁(α, β, γ)
end

function Φ(ℓ::Number, η::Number, ρ::Number)
    return (2*η*ρ)^(ℓ+1) * exp(im*ρ) * M_regularized(ℓ+1+im*η, 2*ℓ+2, -2*im*ρ)
end

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

function F(ℓ::Number, η::Number, ρ::AbstractArray{<:Real})
    return [F(ℓ, η, ρᵢ) for ρᵢ in ρ]
end


export η, C, θ, F, D⁺, D⁻, H⁺, H⁻, F_imag, G, M_regularized, Φ, w
