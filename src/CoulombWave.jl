using SpecialFunctions
using HypergeometricFunctions

"""
    η(a::Number, k::Number)
    η(ϵ::Number)

Coulomb (Sommerfeld) parameter, the second argument of the Coulomb wave
functions such as [`F`](@ref).

- `η(a, k)` returns `1/(a*k)`, where `k` is the wave number and `a` is the
  Bohr radius of the two-body system. `a > 0` describes a repulsive and
  `a < 0` an attractive interaction.
- `η(ϵ)` returns `1/sqrt(ϵ)`, where `ϵ = (a*k)^2` is the energy in units of
  the corresponding Rydberg energy. A negative real `ϵ` throws a
  `DomainError`; pass `complex(ϵ)` for bound-state energies.

Zero arguments throw an `ArgumentError`.
"""
function η(a::Number, k::Number)
    iszero(a) && throw(ArgumentError("a must be nonzero"))
    iszero(k) && throw(ArgumentError("k must be nonzero"))
    return 1 / (a * k)
end

function η(ϵ::Number)
    iszero(ϵ) && throw(ArgumentError("ϵ must be nonzero"))
    return 1 / sqrt(ϵ)
end

"""
    C(ℓ::Number, η::Number)

Coulomb normalization constant (Gamow factor)

    C_ℓ(η) = 2^ℓ * exp(-π*η/2) * sqrt(Γ(ℓ+1+iη) * Γ(ℓ+1-iη)) / Γ(2ℓ+2),

which fixes the behavior `F(ℓ, η, ρ) ≈ C(ℓ, η) * ρ^(ℓ+1)` as `ρ → 0`.
`ℓ` is the angular momentum and `η` the Coulomb parameter (see [`η`](@ref)).
Real arguments return a real value; complex arguments return a complex value.

For complex parameters, the square root of the gamma product is defined by
`exp((loggamma(ℓ + 1 + im * η) + loggamma(ℓ + 1 - im * η)) / 2)`.
This selects a local analytic branch using the principal log-gamma branches,
away from their cuts and poles, and agrees with the positive normalization
for real `ℓ ≥ 0` and real `η`. The full normalization is evaluated in logarithmic
form to avoid overflow or underflow of its individual factors.

See [Eq. (8) of the implementation paper](https://arxiv.org/html/1804.10976v3#S2.E8)
and [DLMF 33.2.5](https://dlmf.nist.gov/33.2.E5).
"""
function C(ℓ::Number, η::Number)
    ℓf, ηf = promote(float(ℓ), float(η))
    logg = (loggamma(ℓf + 1 + im * ηf) + loggamma(ℓf + 1 - im * ηf)) / 2
    value = exp(ℓf * log(2 * one(real(ℓf))) - π * ηf / 2 + logg - loggamma(complex(2 * ℓf + 2)))
    return ℓ isa Real && η isa Real ? real(value) : value
end

"""
    θ(ℓ::Number, η::Number, ρ::Number)

Phase of the Coulomb wave functions at large `ρ`,

    θ_ℓ(η, ρ) = ρ - ℓ*π/2 - η*log(2ρ) + σ_ℓ(η),

where `σ_ℓ(η) = imag(loggamma(ℓ + 1 + im*η))` is the Coulomb phase shift.
As `ρ → ∞` with real arguments, `F ≈ sin(θ)`, `G ≈ cos(θ)`, and
`H⁺ ≈ exp(im*θ)`, `H⁻ ≈ exp(-im*θ)`. See [DLMF 33.2.9](https://dlmf.nist.gov/33.2.E9).
"""
function θ(ℓ::Number, η::Number, ρ::Number)
    logg = loggamma(ℓ + 1 + im * η)
    return ρ - ℓ * π / 2 - η * log(2 * ρ) + imag(logg)
end

"""
    F(ℓ::Number, η::Number, ρ::Number)

Regular Coulomb wave function

    F_ℓ(η, ρ) = C_ℓ(η) * ρ^(ℓ+1) * exp(-iρ) * ₁F₁(ℓ+1-iη; 2ℓ+2; 2iρ),

the solution of the Coulomb wave equation
`u'' + (1 - 2η/ρ - ℓ(ℓ+1)/ρ²) u = 0` that vanishes at `ρ = 0` for `ℓ ≥ 0`
([DLMF 33.2.3](https://dlmf.nist.gov/33.2.E3)).

# Arguments
- `ℓ`: angular momentum. Non-integer and complex values are accepted.
- `η`: Coulomb parameter (see [`η`](@ref)); positive for a repulsive and
  negative for an attractive interaction.
- `ρ`: dimensionless radius `k*r`.

Three real arguments return a real value. If any argument is complex, the
result is complex. The normalization is [`C`](@ref) and the confluent
hypergeometric function is evaluated by HypergeometricFunctions.jl, which
determines the accuracy.

References:
- [Coulomb wave function](https://en.wikipedia.org/wiki/Coulomb_wave_function)
- [Implementation paper](https://arxiv.org/abs/1804.10976)
"""
function F(ℓ::Number, η::Number, ρ::Number)
    ℓc = complex(float(ℓ))
    ηc = complex(float(η))
    ρc = complex(float(ρ))
    return C(ℓc, ηc) * ρc^(ℓc + 1) * exp(-im * ρc) * _₁F₁(complex(ℓc + 1 - im * ηc), complex(2 * ℓc + 2), complex(2 * im * ρc))
end

function F(ℓ::Real, η::Real, ρ::Real)
    return real(F(complex(ℓ), complex(η), complex(ρ)))
end

"""
    D⁺(ℓ::Number, η::Number)

Normalization factor of the outgoing Coulomb wave function [`H⁺`](@ref),

    D⁺_ℓ(η) = (-2i)^(2ℓ+1) * Γ(ℓ+1+iη) / (C_ℓ(η) * Γ(2ℓ+2)),

with `C_ℓ(η)` given by [`C`](@ref). The result is complex.
"""
function D⁺(ℓ::Number, η::Number)
    return (-2 * im)^(2 * ℓ + 1) * gamma(ℓ + 1 + im * η) / (C(ℓ, η) * gamma(2 * ℓ + 2))
end

"""
    D⁻(ℓ::Number, η::Number)

Normalization factor of the incoming Coulomb wave function [`H⁻`](@ref),

    D⁻_ℓ(η) = (2i)^(2ℓ+1) * Γ(ℓ+1-iη) / (C_ℓ(η) * Γ(2ℓ+2)),

with `C_ℓ(η)` given by [`C`](@ref). The result is complex; for real `ℓ` and
`η` it is the complex conjugate of [`D⁺`](@ref).
"""
function D⁻(ℓ::Number, η::Number)
    return (2 * im)^(2 * ℓ + 1) * gamma(ℓ + 1 - im * η) / (C(ℓ, η) * gamma(2 * ℓ + 2))
end

"""
    H⁺(ℓ::Number, η::Number, ρ::Number)

Outgoing Coulomb wave function

    H⁺_ℓ(η, ρ) = D⁺_ℓ(η) * ρ^(ℓ+1) * exp(iρ) * U(ℓ+1+iη, 2ℓ+2, -2iρ),

where `U` is Tricomi's confluent hypergeometric function and `D⁺` is the
factor [`D⁺`](@ref). The result is complex. For real arguments
`H⁺ = G + im*F`, and `H⁺ ≈ exp(im*θ)` as `ρ → ∞` with the phase [`θ`](@ref).
The arguments have the same meaning as in [`F`](@ref); `ρ` must be nonzero.

References:
- [Coulomb wave function](https://en.wikipedia.org/wiki/Coulomb_wave_function)
- [Implementation paper](https://arxiv.org/abs/1804.10976)
"""
function H⁺(ℓ::Number, η::Number, ρ::Number)
    return D⁺(ℓ, η) * ρ^(ℓ + 1) * exp(+im * ρ) * HypergeometricFunctions.U(ℓ + 1 + im * η, 2 * ℓ + 2, -2 * im * ρ)
end

"""
    H⁻(ℓ::Number, η::Number, ρ::Number)

Incoming Coulomb wave function

    H⁻_ℓ(η, ρ) = D⁻_ℓ(η) * ρ^(ℓ+1) * exp(-iρ) * U(ℓ+1-iη, 2ℓ+2, 2iρ),

where `U` is Tricomi's confluent hypergeometric function and `D⁻` is the
factor [`D⁻`](@ref). The result is complex. For real arguments
`H⁻ = G - im*F`, and `H⁻ ≈ exp(-im*θ)` as `ρ → ∞` with the phase [`θ`](@ref).
The arguments have the same meaning as in [`F`](@ref); `ρ` must be nonzero.

References:
- [Coulomb wave function](https://en.wikipedia.org/wiki/Coulomb_wave_function)
- [Implementation paper](https://arxiv.org/abs/1804.10976)
"""
function H⁻(ℓ::Number, η::Number, ρ::Number)
    return D⁻(ℓ, η) * ρ^(ℓ + 1) * exp(-im * ρ) * HypergeometricFunctions.U(ℓ + 1 - im * η, 2 * ℓ + 2, +2 * im * ρ)
end

"""
    F_imag(ℓ::Number, η::Number, ρ::Number)

Regular Coulomb wave function computed from the outgoing and incoming
solutions, `(H⁺(ℓ, η, ρ) - H⁻(ℓ, η, ρ)) / (2im)`.

For real arguments this is the imaginary part of [`H⁺`](@ref) and agrees
with [`F`](@ref) up to rounding, but the result is returned as a complex
number and, unlike `F`, the evaluation requires `ρ ≠ 0`.
"""
function F_imag(ℓ::Number, η::Number, ρ::Number)
    return (H⁺(ℓ, η, ρ) - H⁻(ℓ, η, ρ)) / (2 * im)
end

"""
    G(ℓ::Number, η::Number, ρ::Number)

Irregular Coulomb wave function

    G_ℓ(η, ρ) = (H⁺_ℓ(η, ρ) + H⁻_ℓ(η, ρ)) / 2,

the solution of the Coulomb wave equation that is linearly independent of
[`F`](@ref), with Wronskian `F'G - FG' = 1` and `G ≈ cos(θ)` as `ρ → ∞`
(see [`θ`](@ref)). It is singular at `ρ = 0`.

The arguments have the same meaning as in [`F`](@ref). Three real arguments
return a real value; otherwise the result is complex.

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

Regularized confluent hypergeometric function
`₁F₁(α; β; γ) / Γ(β)`, with parameters `α` and `β` and argument `γ`.

The value is computed as the quotient of the two factors, so it is not
finite at the poles of `Γ(β)`, `β = 0, -1, -2, …`, where the mathematical
function has a finite limit.
"""
function M_regularized(α::Number, β::Number, γ::Number)
    return 1 / gamma(β) * HypergeometricFunctions._₁F₁(α, β, γ)
end

"""
    Φ(ℓ::Number, η::Number, ρ::Number)

Modified regular Coulomb function

    Φ_ℓ(η, ρ) = (2ηρ)^(ℓ+1) * exp(iρ) * M_regularized(ℓ+1+iη, 2ℓ+2, -2iρ),

which equals `F(ℓ, η, ρ)` up to the `ρ`-independent factor
`(2η)^(ℓ+1) / (C(ℓ, η) * Γ(2ℓ+2))`. See [`M_regularized`](@ref) and the
[implementation paper](https://arxiv.org/abs/1804.10976). The arguments have
the same meaning as in [`F`](@ref). The result is complex.
"""
function Φ(ℓ::Number, η::Number, ρ::Number)
    return (2 * η * ρ)^(ℓ + 1) * exp(im * ρ) * M_regularized(ℓ + 1 + im * η, 2 * ℓ + 2, -2 * im * ρ)
end

"""
    w(ℓ::Integer, η::Number)
    w(ℓ::Number, η::Number)

Product

    w_ℓ(η) = ∏ⱼ (1 + j²/η²),

taken over `j = 0, 1, …, ℓ` for an integer `ℓ` and over
`j = 1/2, 3/2, …, ℓ` for a half-integer `ℓ`. It appears in the connection
formulas of the [implementation paper](https://arxiv.org/abs/1804.10976)
(see [`Ψ`](@ref)).

An integer order must be passed as an `Integer`: `w(1, η)` is accepted,
while `w(1.0, η)` and any other non-half-integer value throw an
`ArgumentError`.
"""
function w(ℓ::Integer, η::Number)
    result = one(η)
    for j in 0:ℓ
        result *= 1 + j^2 / η^2
    end
    return result
end

function w(ℓ::Number, η::Number)
    T = typeof(float(ℓ))
    if isapprox(mod(ℓ - T(0.5), one(T)), zero(T); atol = eps(T) * 100)
        result = one(η)
        j = T(0.5)
        while j <= ℓ
            result *= 1 + j^2 / η^2
            j += 1
        end
        return result
    else
        throw(ArgumentError("ℓ must be either an integer or half-integer (1/2, 3/2, 5/2, ...)"))
    end
end

"""
    w_plus(ℓ::Number, η::Number)

Gamma-function ratio

    w⁺_ℓ(η) = Γ(ℓ+1+iη) / ((iη)^(2ℓ+1) * Γ(-ℓ+iη))

used in the connection formulas of the
[implementation paper](https://arxiv.org/abs/1804.10976). The result is
complex; `η` must be nonzero.
"""
function w_plus(ℓ::Number, η::Number)
    return gamma(ℓ + 1 + im * η) / ((im * η)^(2 * ℓ + 1) * gamma(-ℓ + im * η))
end

"""
    w_minus(ℓ::Number, η::Number)

Gamma-function ratio

    w⁻_ℓ(η) = Γ(ℓ+1-iη) / ((-iη)^(2ℓ+1) * Γ(-ℓ-iη)),

the counterpart of [`w_plus`](@ref) with `i → -i`. The result is complex;
`η` must be nonzero.
"""
function w_minus(ℓ::Number, η::Number)
    return gamma(ℓ + 1 - im * η) / ((-im * η)^(2 * ℓ + 1) * gamma(-ℓ - im * η))
end

"""
    h_plus(ℓ::Number, η::Number)

Digamma combination

    h⁺_ℓ(η) = (ψ(ℓ+1+iη) + ψ(-ℓ+iη)) / 2 - log(iη),

where `ψ` is the digamma function and `log` is the principal branch. It is
used in the connection formulas of the
[implementation paper](https://arxiv.org/abs/1804.10976). The result is
complex; `η` must be nonzero.
"""
function h_plus(ℓ::Number, η::Number)
    return (digamma(ℓ + 1 + im * η) + digamma(-ℓ + im * η)) / 2 - log(im * η)
end

"""
    h_minus(ℓ::Number, η::Number)

Digamma combination

    h⁻_ℓ(η) = (ψ(ℓ+1-iη) + ψ(-ℓ-iη)) / 2 - log(-iη),

the counterpart of [`h_plus`](@ref) with `i → -i`. The result is complex;
`η` must be nonzero.
"""
function h_minus(ℓ::Number, η::Number)
    return (digamma(ℓ + 1 - im * η) + digamma(-ℓ - im * η)) / 2 - log(-im * η)
end

"""
    g(ℓ::Number, η::Number)

Real part of the digamma combination

    g_ℓ(η) = real((ψ(ℓ+1+iη) + ψ(ℓ+1-iη)) / 2) - log(|η|),

where `ψ` is the digamma function. It is used in the connection formulas of
the [implementation paper](https://arxiv.org/abs/1804.10976). `η` must be
nonzero.
"""
function g(ℓ::Number, η::Number)
    x = (digamma(ℓ + 1 + im * η) + digamma(ℓ + 1 - im * η)) / 2 - log(abs(η))
    return real(x)
end

"""
    Φ_dot(ℓ::Number, η::Number, ρ::Number; h=nothing)

Derivative `∂Φ_ℓ(η, ρ)/∂ℓ` of [`Φ`](@ref) with respect to the angular
momentum, approximated by the central difference
`(Φ(ℓ + h, η, ρ) - Φ(ℓ - h, η, ρ)) / (2h)`.

`h` is the step size; it defaults to `cbrt(eps(T))`, where `T` is the
floating-point type of `ℓ`. With that step the truncation and rounding errors
are both of order `eps(T)^(2/3)` relative to the scale of `Φ` (about `4e-11`
for `Float64`), so the result has fewer correct digits than `Φ` itself.
"""
function Φ_dot(ℓ::Number, η::Number, ρ::Number; h = nothing)
    h_eff = h === nothing ? cbrt(eps(real(float(typeof(ℓ))))) : h
    return (Φ(ℓ + h_eff, η, ρ) - Φ(ℓ - h_eff, η, ρ)) / (2 * h_eff)
end

"""
    F_dot(ℓ::Number, η::Number, ρ::Number; h=nothing)

Derivative `∂F_ℓ(η, ρ)/∂ℓ` of [`F`](@ref) with respect to the angular
momentum, approximated by the central difference
`(F(ℓ + h, η, ρ) - F(ℓ - h, η, ρ)) / (2h)`.

`h` is the step size; it defaults to `cbrt(eps(T))`, where `T` is the
floating-point type of `ℓ`. The accuracy is that of [`Φ_dot`](@ref).
"""
function F_dot(ℓ::Number, η::Number, ρ::Number; h = nothing)
    h_eff = h === nothing ? cbrt(eps(real(float(typeof(ℓ))))) : h
    return (F(ℓ + h_eff, η, ρ) - F(ℓ - h_eff, η, ρ)) / (2 * h_eff)
end

"""
    Ψ(ℓ::Number, η::Number, ρ::Number; h=nothing)

Combination of angular-momentum derivatives of [`Φ`](@ref),

    Ψ_ℓ(η, ρ) = (w_ℓ(η) * ∂Φ_ℓ/∂ℓ + ∂Φ_{-ℓ-1}/∂ℓ) / 2,

with `w_ℓ(η)` given by [`w`](@ref) and both derivatives computed by
[`Φ_dot`](@ref). See the
[implementation paper](https://arxiv.org/abs/1804.10976) for its role in the
connection formulas.

`ℓ` must be an `Integer` or a half-integer, as required by `w`. `h` is the
finite-difference step passed to `Φ_dot`, whose accuracy limits apply.
"""
function Ψ(ℓ::Number, η::Number, ρ::Number; h = nothing)
    h_eff = h === nothing ? cbrt(eps(real(float(typeof(ℓ))))) : h
    return w(ℓ, η) * Φ_dot(ℓ, η, ρ; h = h_eff) / 2 + Φ_dot(-ℓ - 1, η, ρ; h = h_eff) / 2
end

"""
    I(ℓ::Number, η::Number, ρ::Number; h=nothing)

Rescaled form of [`Ψ`](@ref),

    I_ℓ(η, ρ) = C_ℓ(η) * Γ(2ℓ+2) / (2η)^(ℓ+1) * Ψ_ℓ(η, ρ),

which applies to `Ψ` the inverse of the factor that relates [`Φ`](@ref) to
[`F`](@ref). `C_ℓ(η)` is [`C`](@ref).

`ℓ` must be an `Integer` or a half-integer, and `η` must be nonzero. `h` is
the finite-difference step passed to [`Φ_dot`](@ref), whose accuracy limits
apply. This function shares its name with `LinearAlgebra.I`; qualify it as
`FewSpecialFunctions.I` when both packages are loaded.
"""
function I(ℓ::Number, η::Number, ρ::Number; h = nothing)
    h_eff = h === nothing ? cbrt(eps(real(float(typeof(ℓ))))) : h
    return C(ℓ, η) * gamma(2 * ℓ + 2) / ((2 * η)^(ℓ + 1)) * Ψ(ℓ, η, ρ; h = h_eff)
end


export η, C, θ, F, D⁺, D⁻, H⁺, H⁻, F_imag, G, M_regularized, Φ, w, w_plus, w_minus, h_plus, h_minus, g, Φ_dot, F_dot, Ψ, I
