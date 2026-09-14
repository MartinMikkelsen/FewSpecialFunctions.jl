module FewSpecialFunctionsForwardDiffExt

using FewSpecialFunctions
using ForwardDiff
using ForwardDiff: Dual, partials, value

import FewSpecialFunctions:
    η, C, θ, F, D⁺, D⁻, H⁺, H⁻, F_imag, G, M_regularized, Φ, w,
    debye_function,
    fresnel, FresnelS, FresnelC, FresnelE,
    dawson,
    voigt,
    Clausen,
    FermiDiracIntegral, FermiDiracIntegralNorm,
    BoseEinsteinIntegral, BoseEinsteinIntegralNorm,
    MarcumQ, dQdb,
    U, V, W, dU, dV, dW, U_scaled, V_scaled,
    WhittakerM, WhittakerW, dWhittakerM, dWhittakerW

_fd_step(x::T) where {T <: Real} = cbrt(eps(T)) * (abs(x) + one(T))
_fd_deriv(f, x) = (h = _fd_step(x); (f(x + h) - f(x - h)) / (2h))

# For complex-valued functions: split the derivative into real and imaginary Dual parts.
_complex_dual(::Type{T}, y, dy, p) where {T} =
    complex(Dual{T}(real(y), real(dy) * p), Dual{T}(imag(y), imag(dy) * p))

# ── Coulomb: η ─────────────────────────────────────────────────────────────────
# dη/dϵ = -1 / (2 ϵ^(3/2))

function η(ϵ::Dual{T}) where {T}
    ϵv = value(ϵ)
    return Dual{T}(η(ϵv), -one(ϵv) / (2 * ϵv * sqrt(ϵv)) * partials(ϵ))
end

# dη/da = -1/(a²k),  dη/dk = -1/(ak²)

function η(a::Dual{T}, k::Real) where {T}
    av = value(a)
    return Dual{T}(η(av, k), -one(av) / (av^2 * k) * partials(a))
end

function η(a::Real, k::Dual{T}) where {T}
    kv = value(k)
    return Dual{T}(η(a, kv), -one(kv) / (a * kv^2) * partials(k))
end

function η(a::Dual{T}, k::Dual{T}) where {T}
    av, kv = value(a), value(k)
    return Dual{T}(η(av, kv), (-one(av) / (av^2 * kv)) * partials(a) + (-one(kv) / (av * kv^2)) * partials(k))
end

# ── Coulomb: C (real-valued, FD) ───────────────────────────────────────────────

function C(ℓ::Dual{T}, η_::Real) where {T}
    ℓv = value(ℓ)
    return Dual{T}(C(ℓv, η_), _fd_deriv(x -> C(x, η_), ℓv) * partials(ℓ))
end

function C(ℓ::Real, η_::Dual{T}) where {T}
    ηv = value(η_)
    return Dual{T}(C(ℓ, ηv), _fd_deriv(x -> C(ℓ, x), ηv) * partials(η_))
end

function C(ℓ::Dual{T}, η_::Dual{T}) where {T}
    ℓv, ηv = value(ℓ), value(η_)
    y = C(ℓv, ηv)
    return Dual{T}(y, _fd_deriv(x -> C(x, ηv), ℓv) * partials(ℓ) + _fd_deriv(x -> C(ℓv, x), ηv) * partials(η_))
end

# ── Voigt (real-valued, FD) ───────────────────────────────────────────────────

function voigt(x::Dual{T}, y::Real) where {T}
    xv = value(x)
    return Dual{T}(voigt(xv, y), _fd_deriv(t -> voigt(t, y), xv) * partials(x))
end

function voigt(x::Real, y::Dual{T}) where {T}
    yv = value(y)
    return Dual{T}(voigt(x, yv), _fd_deriv(t -> voigt(x, t), yv) * partials(y))
end

function voigt(x::Dual{T}, y::Dual{T}) where {T}
    xv, yv = value(x), value(y)
    z = voigt(xv, yv)
    return Dual{T}(z, _fd_deriv(t -> voigt(t, yv), xv) * partials(x) + _fd_deriv(t -> voigt(xv, t), yv) * partials(y))
end

# ── Coulomb: D⁺, D⁻ (complex-valued, FD) ──────────────────────────────────────

function D⁺(ℓ::Dual{T}, η_::Real) where {T}
    ℓv = value(ℓ)
    y = D⁺(ℓv, η_)
    return _complex_dual(T, y, _fd_deriv(x -> D⁺(x, η_), ℓv), partials(ℓ))
end

function D⁺(ℓ::Real, η_::Dual{T}) where {T}
    ηv = value(η_)
    y = D⁺(ℓ, ηv)
    return _complex_dual(T, y, _fd_deriv(x -> D⁺(ℓ, x), ηv), partials(η_))
end

function D⁺(ℓ::Dual{T}, η_::Dual{T}) where {T}
    ℓv, ηv = value(ℓ), value(η_)
    y = D⁺(ℓv, ηv)
    dℓ = _fd_deriv(x -> D⁺(x, ηv), ℓv)
    dη = _fd_deriv(x -> D⁺(ℓv, x), ηv)
    return complex(
        Dual{T}(real(y), real(dℓ) * partials(ℓ) + real(dη) * partials(η_)),
        Dual{T}(imag(y), imag(dℓ) * partials(ℓ) + imag(dη) * partials(η_)),
    )
end

function D⁻(ℓ::Dual{T}, η_::Real) where {T}
    ℓv = value(ℓ)
    y = D⁻(ℓv, η_)
    return _complex_dual(T, y, _fd_deriv(x -> D⁻(x, η_), ℓv), partials(ℓ))
end

function D⁻(ℓ::Real, η_::Dual{T}) where {T}
    ηv = value(η_)
    y = D⁻(ℓ, ηv)
    return _complex_dual(T, y, _fd_deriv(x -> D⁻(ℓ, x), ηv), partials(η_))
end

function D⁻(ℓ::Dual{T}, η_::Dual{T}) where {T}
    ℓv, ηv = value(ℓ), value(η_)
    y = D⁻(ℓv, ηv)
    dℓ = _fd_deriv(x -> D⁻(x, ηv), ℓv)
    dη = _fd_deriv(x -> D⁻(ℓv, x), ηv)
    return complex(
        Dual{T}(real(y), real(dℓ) * partials(ℓ) + real(dη) * partials(η_)),
        Dual{T}(imag(y), imag(dℓ) * partials(ℓ) + imag(dη) * partials(η_)),
    )
end

# ── Coulomb: θ (real-valued; analytic dθ/dρ = 1 - η/ρ, FD for ℓ and η) ───────

function θ(ℓ::Dual{T}, η_::Real, ρ::Real) where {T}
    ℓv = value(ℓ)
    return Dual{T}(θ(ℓv, η_, ρ), _fd_deriv(x -> θ(x, η_, ρ), ℓv) * partials(ℓ))
end

function θ(ℓ::Real, η_::Dual{T}, ρ::Real) where {T}
    ηv = value(η_)
    return Dual{T}(θ(ℓ, ηv, ρ), _fd_deriv(x -> θ(ℓ, x, ρ), ηv) * partials(η_))
end

function θ(ℓ::Real, η_::Real, ρ::Dual{T}) where {T}
    ρv = value(ρ)
    return Dual{T}(θ(ℓ, η_, ρv), (one(ρv) - η_ / ρv) * partials(ρ))
end

function θ(ℓ::Dual{T}, η_::Dual{T}, ρ::Real) where {T}
    ℓv, ηv = value(ℓ), value(η_)
    y = θ(ℓv, ηv, ρ)
    return Dual{T}(y, _fd_deriv(x -> θ(x, ηv, ρ), ℓv) * partials(ℓ) + _fd_deriv(x -> θ(ℓv, x, ρ), ηv) * partials(η_))
end

function θ(ℓ::Dual{T}, η_::Real, ρ::Dual{T}) where {T}
    ℓv, ρv = value(ℓ), value(ρ)
    y = θ(ℓv, η_, ρv)
    return Dual{T}(y, _fd_deriv(x -> θ(x, η_, ρv), ℓv) * partials(ℓ) + (one(ρv) - η_ / ρv) * partials(ρ))
end

function θ(ℓ::Real, η_::Dual{T}, ρ::Dual{T}) where {T}
    ηv, ρv = value(η_), value(ρ)
    y = θ(ℓ, ηv, ρv)
    return Dual{T}(y, _fd_deriv(x -> θ(ℓ, x, ρv), ηv) * partials(η_) + (one(ρv) - ηv / ρv) * partials(ρ))
end

function θ(ℓ::Dual{T}, η_::Dual{T}, ρ::Dual{T}) where {T}
    ℓv, ηv, ρv = value(ℓ), value(η_), value(ρ)
    y = θ(ℓv, ηv, ρv)
    return Dual{T}(
        y,
        _fd_deriv(x -> θ(x, ηv, ρv), ℓv) * partials(ℓ) +
            _fd_deriv(x -> θ(ℓv, x, ρv), ηv) * partials(η_) +
            (one(ρv) - ηv / ρv) * partials(ρ),
    )
end

# ── Coulomb: F (real-valued, FD) ───────────────────────────────────────────────

function F(ℓ::Dual{T}, η_::Real, ρ::Real) where {T}
    ℓv = value(ℓ)
    return Dual{T}(F(ℓv, η_, ρ), _fd_deriv(x -> F(x, η_, ρ), ℓv) * partials(ℓ))
end

function F(ℓ::Real, η_::Dual{T}, ρ::Real) where {T}
    ηv = value(η_)
    return Dual{T}(F(ℓ, ηv, ρ), _fd_deriv(x -> F(ℓ, x, ρ), ηv) * partials(η_))
end

function F(ℓ::Real, η_::Real, ρ::Dual{T}) where {T}
    ρv = value(ρ)
    return Dual{T}(F(ℓ, η_, ρv), _fd_deriv(x -> F(ℓ, η_, x), ρv) * partials(ρ))
end

function F(ℓ::Dual{T}, η_::Dual{T}, ρ::Real) where {T}
    ℓv, ηv = value(ℓ), value(η_)
    y = F(ℓv, ηv, ρ)
    return Dual{T}(y, _fd_deriv(x -> F(x, ηv, ρ), ℓv) * partials(ℓ) + _fd_deriv(x -> F(ℓv, x, ρ), ηv) * partials(η_))
end

function F(ℓ::Dual{T}, η_::Real, ρ::Dual{T}) where {T}
    ℓv, ρv = value(ℓ), value(ρ)
    y = F(ℓv, η_, ρv)
    return Dual{T}(y, _fd_deriv(x -> F(x, η_, ρv), ℓv) * partials(ℓ) + _fd_deriv(x -> F(ℓv, η_, x), ρv) * partials(ρ))
end

function F(ℓ::Real, η_::Dual{T}, ρ::Dual{T}) where {T}
    ηv, ρv = value(η_), value(ρ)
    y = F(ℓ, ηv, ρv)
    return Dual{T}(y, _fd_deriv(x -> F(ℓ, x, ρv), ηv) * partials(η_) + _fd_deriv(x -> F(ℓ, ηv, x), ρv) * partials(ρ))
end

function F(ℓ::Dual{T}, η_::Dual{T}, ρ::Dual{T}) where {T}
    ℓv, ηv, ρv = value(ℓ), value(η_), value(ρ)
    y = F(ℓv, ηv, ρv)
    return Dual{T}(
        y,
        _fd_deriv(x -> F(x, ηv, ρv), ℓv) * partials(ℓ) +
            _fd_deriv(x -> F(ℓv, x, ρv), ηv) * partials(η_) +
            _fd_deriv(x -> F(ℓv, ηv, x), ρv) * partials(ρ),
    )
end

# ── Coulomb: H⁺ (complex-valued, FD) ──────────────────────────────────────────

function H⁺(ℓ::Dual{T}, η_::Real, ρ::Real) where {T}
    ℓv = value(ℓ)
    y = H⁺(ℓv, η_, ρ)
    return _complex_dual(T, y, _fd_deriv(x -> H⁺(x, η_, ρ), ℓv), partials(ℓ))
end

function H⁺(ℓ::Real, η_::Dual{T}, ρ::Real) where {T}
    ηv = value(η_)
    y = H⁺(ℓ, ηv, ρ)
    return _complex_dual(T, y, _fd_deriv(x -> H⁺(ℓ, x, ρ), ηv), partials(η_))
end

function H⁺(ℓ::Real, η_::Real, ρ::Dual{T}) where {T}
    ρv = value(ρ)
    y = H⁺(ℓ, η_, ρv)
    return _complex_dual(T, y, _fd_deriv(x -> H⁺(ℓ, η_, x), ρv), partials(ρ))
end

function H⁺(ℓ::Dual{T}, η_::Dual{T}, ρ::Real) where {T}
    ℓv, ηv = value(ℓ), value(η_)
    y = H⁺(ℓv, ηv, ρ)
    dℓ = _fd_deriv(x -> H⁺(x, ηv, ρ), ℓv)
    dη = _fd_deriv(x -> H⁺(ℓv, x, ρ), ηv)
    return complex(Dual{T}(real(y), real(dℓ) * partials(ℓ) + real(dη) * partials(η_)), Dual{T}(imag(y), imag(dℓ) * partials(ℓ) + imag(dη) * partials(η_)))
end

function H⁺(ℓ::Dual{T}, η_::Real, ρ::Dual{T}) where {T}
    ℓv, ρv = value(ℓ), value(ρ)
    y = H⁺(ℓv, η_, ρv)
    dℓ = _fd_deriv(x -> H⁺(x, η_, ρv), ℓv)
    dρ = _fd_deriv(x -> H⁺(ℓv, η_, x), ρv)
    return complex(Dual{T}(real(y), real(dℓ) * partials(ℓ) + real(dρ) * partials(ρ)), Dual{T}(imag(y), imag(dℓ) * partials(ℓ) + imag(dρ) * partials(ρ)))
end

function H⁺(ℓ::Real, η_::Dual{T}, ρ::Dual{T}) where {T}
    ηv, ρv = value(η_), value(ρ)
    y = H⁺(ℓ, ηv, ρv)
    dη = _fd_deriv(x -> H⁺(ℓ, x, ρv), ηv)
    dρ = _fd_deriv(x -> H⁺(ℓ, ηv, x), ρv)
    return complex(Dual{T}(real(y), real(dη) * partials(η_) + real(dρ) * partials(ρ)), Dual{T}(imag(y), imag(dη) * partials(η_) + imag(dρ) * partials(ρ)))
end

function H⁺(ℓ::Dual{T}, η_::Dual{T}, ρ::Dual{T}) where {T}
    ℓv, ηv, ρv = value(ℓ), value(η_), value(ρ)
    y = H⁺(ℓv, ηv, ρv)
    dℓ = _fd_deriv(x -> H⁺(x, ηv, ρv), ℓv)
    dη = _fd_deriv(x -> H⁺(ℓv, x, ρv), ηv)
    dρ = _fd_deriv(x -> H⁺(ℓv, ηv, x), ρv)
    return complex(
        Dual{T}(real(y), real(dℓ) * partials(ℓ) + real(dη) * partials(η_) + real(dρ) * partials(ρ)),
        Dual{T}(imag(y), imag(dℓ) * partials(ℓ) + imag(dη) * partials(η_) + imag(dρ) * partials(ρ)),
    )
end

# ── Coulomb: H⁻ (complex-valued, FD) ──────────────────────────────────────────

function H⁻(ℓ::Dual{T}, η_::Real, ρ::Real) where {T}
    ℓv = value(ℓ)
    y = H⁻(ℓv, η_, ρ)
    return _complex_dual(T, y, _fd_deriv(x -> H⁻(x, η_, ρ), ℓv), partials(ℓ))
end

function H⁻(ℓ::Real, η_::Dual{T}, ρ::Real) where {T}
    ηv = value(η_)
    y = H⁻(ℓ, ηv, ρ)
    return _complex_dual(T, y, _fd_deriv(x -> H⁻(ℓ, x, ρ), ηv), partials(η_))
end

function H⁻(ℓ::Real, η_::Real, ρ::Dual{T}) where {T}
    ρv = value(ρ)
    y = H⁻(ℓ, η_, ρv)
    return _complex_dual(T, y, _fd_deriv(x -> H⁻(ℓ, η_, x), ρv), partials(ρ))
end

function H⁻(ℓ::Dual{T}, η_::Dual{T}, ρ::Real) where {T}
    ℓv, ηv = value(ℓ), value(η_)
    y = H⁻(ℓv, ηv, ρ)
    dℓ = _fd_deriv(x -> H⁻(x, ηv, ρ), ℓv)
    dη = _fd_deriv(x -> H⁻(ℓv, x, ρ), ηv)
    return complex(Dual{T}(real(y), real(dℓ) * partials(ℓ) + real(dη) * partials(η_)), Dual{T}(imag(y), imag(dℓ) * partials(ℓ) + imag(dη) * partials(η_)))
end

function H⁻(ℓ::Dual{T}, η_::Real, ρ::Dual{T}) where {T}
    ℓv, ρv = value(ℓ), value(ρ)
    y = H⁻(ℓv, η_, ρv)
    dℓ = _fd_deriv(x -> H⁻(x, η_, ρv), ℓv)
    dρ = _fd_deriv(x -> H⁻(ℓv, η_, x), ρv)
    return complex(Dual{T}(real(y), real(dℓ) * partials(ℓ) + real(dρ) * partials(ρ)), Dual{T}(imag(y), imag(dℓ) * partials(ℓ) + imag(dρ) * partials(ρ)))
end

function H⁻(ℓ::Real, η_::Dual{T}, ρ::Dual{T}) where {T}
    ηv, ρv = value(η_), value(ρ)
    y = H⁻(ℓ, ηv, ρv)
    dη = _fd_deriv(x -> H⁻(ℓ, x, ρv), ηv)
    dρ = _fd_deriv(x -> H⁻(ℓ, ηv, x), ρv)
    return complex(Dual{T}(real(y), real(dη) * partials(η_) + real(dρ) * partials(ρ)), Dual{T}(imag(y), imag(dη) * partials(η_) + imag(dρ) * partials(ρ)))
end

function H⁻(ℓ::Dual{T}, η_::Dual{T}, ρ::Dual{T}) where {T}
    ℓv, ηv, ρv = value(ℓ), value(η_), value(ρ)
    y = H⁻(ℓv, ηv, ρv)
    dℓ = _fd_deriv(x -> H⁻(x, ηv, ρv), ℓv)
    dη = _fd_deriv(x -> H⁻(ℓv, x, ρv), ηv)
    dρ = _fd_deriv(x -> H⁻(ℓv, ηv, x), ρv)
    return complex(
        Dual{T}(real(y), real(dℓ) * partials(ℓ) + real(dη) * partials(η_) + real(dρ) * partials(ρ)),
        Dual{T}(imag(y), imag(dℓ) * partials(ℓ) + imag(dη) * partials(η_) + imag(dρ) * partials(ρ)),
    )
end

# ── Coulomb: F_imag (complex-valued, FD) ───────────────────────────────────────

function F_imag(ℓ::Dual{T}, η_::Real, ρ::Real) where {T}
    ℓv = value(ℓ)
    y = F_imag(ℓv, η_, ρ)
    return _complex_dual(T, y, _fd_deriv(x -> F_imag(x, η_, ρ), ℓv), partials(ℓ))
end

function F_imag(ℓ::Real, η_::Dual{T}, ρ::Real) where {T}
    ηv = value(η_)
    y = F_imag(ℓ, ηv, ρ)
    return _complex_dual(T, y, _fd_deriv(x -> F_imag(ℓ, x, ρ), ηv), partials(η_))
end

function F_imag(ℓ::Real, η_::Real, ρ::Dual{T}) where {T}
    ρv = value(ρ)
    y = F_imag(ℓ, η_, ρv)
    return _complex_dual(T, y, _fd_deriv(x -> F_imag(ℓ, η_, x), ρv), partials(ρ))
end

function F_imag(ℓ::Dual{T}, η_::Dual{T}, ρ::Real) where {T}
    ℓv, ηv = value(ℓ), value(η_)
    y = F_imag(ℓv, ηv, ρ)
    dℓ = _fd_deriv(x -> F_imag(x, ηv, ρ), ℓv)
    dη = _fd_deriv(x -> F_imag(ℓv, x, ρ), ηv)
    return complex(Dual{T}(real(y), real(dℓ) * partials(ℓ) + real(dη) * partials(η_)), Dual{T}(imag(y), imag(dℓ) * partials(ℓ) + imag(dη) * partials(η_)))
end

function F_imag(ℓ::Dual{T}, η_::Real, ρ::Dual{T}) where {T}
    ℓv, ρv = value(ℓ), value(ρ)
    y = F_imag(ℓv, η_, ρv)
    dℓ = _fd_deriv(x -> F_imag(x, η_, ρv), ℓv)
    dρ = _fd_deriv(x -> F_imag(ℓv, η_, x), ρv)
    return complex(Dual{T}(real(y), real(dℓ) * partials(ℓ) + real(dρ) * partials(ρ)), Dual{T}(imag(y), imag(dℓ) * partials(ℓ) + imag(dρ) * partials(ρ)))
end

function F_imag(ℓ::Real, η_::Dual{T}, ρ::Dual{T}) where {T}
    ηv, ρv = value(η_), value(ρ)
    y = F_imag(ℓ, ηv, ρv)
    dη = _fd_deriv(x -> F_imag(ℓ, x, ρv), ηv)
    dρ = _fd_deriv(x -> F_imag(ℓ, ηv, x), ρv)
    return complex(Dual{T}(real(y), real(dη) * partials(η_) + real(dρ) * partials(ρ)), Dual{T}(imag(y), imag(dη) * partials(η_) + imag(dρ) * partials(ρ)))
end

function F_imag(ℓ::Dual{T}, η_::Dual{T}, ρ::Dual{T}) where {T}
    ℓv, ηv, ρv = value(ℓ), value(η_), value(ρ)
    y = F_imag(ℓv, ηv, ρv)
    dℓ = _fd_deriv(x -> F_imag(x, ηv, ρv), ℓv)
    dη = _fd_deriv(x -> F_imag(ℓv, x, ρv), ηv)
    dρ = _fd_deriv(x -> F_imag(ℓv, ηv, x), ρv)
    return complex(
        Dual{T}(real(y), real(dℓ) * partials(ℓ) + real(dη) * partials(η_) + real(dρ) * partials(ρ)),
        Dual{T}(imag(y), imag(dℓ) * partials(ℓ) + imag(dη) * partials(η_) + imag(dρ) * partials(ρ)),
    )
end

# ── Coulomb: Φ (complex-valued, FD) ────────────────────────────────────────────

function Φ(ℓ::Dual{T}, η_::Real, ρ::Real) where {T}
    ℓv = value(ℓ)
    y = Φ(ℓv, η_, ρ)
    return _complex_dual(T, y, _fd_deriv(x -> Φ(x, η_, ρ), ℓv), partials(ℓ))
end

function Φ(ℓ::Real, η_::Dual{T}, ρ::Real) where {T}
    ηv = value(η_)
    y = Φ(ℓ, ηv, ρ)
    return _complex_dual(T, y, _fd_deriv(x -> Φ(ℓ, x, ρ), ηv), partials(η_))
end

function Φ(ℓ::Real, η_::Real, ρ::Dual{T}) where {T}
    ρv = value(ρ)
    y = Φ(ℓ, η_, ρv)
    return _complex_dual(T, y, _fd_deriv(x -> Φ(ℓ, η_, x), ρv), partials(ρ))
end

function Φ(ℓ::Dual{T}, η_::Dual{T}, ρ::Real) where {T}
    ℓv, ηv = value(ℓ), value(η_)
    y = Φ(ℓv, ηv, ρ)
    dℓ = _fd_deriv(x -> Φ(x, ηv, ρ), ℓv)
    dη = _fd_deriv(x -> Φ(ℓv, x, ρ), ηv)
    return complex(Dual{T}(real(y), real(dℓ) * partials(ℓ) + real(dη) * partials(η_)), Dual{T}(imag(y), imag(dℓ) * partials(ℓ) + imag(dη) * partials(η_)))
end

function Φ(ℓ::Dual{T}, η_::Real, ρ::Dual{T}) where {T}
    ℓv, ρv = value(ℓ), value(ρ)
    y = Φ(ℓv, η_, ρv)
    dℓ = _fd_deriv(x -> Φ(x, η_, ρv), ℓv)
    dρ = _fd_deriv(x -> Φ(ℓv, η_, x), ρv)
    return complex(Dual{T}(real(y), real(dℓ) * partials(ℓ) + real(dρ) * partials(ρ)), Dual{T}(imag(y), imag(dℓ) * partials(ℓ) + imag(dρ) * partials(ρ)))
end

function Φ(ℓ::Real, η_::Dual{T}, ρ::Dual{T}) where {T}
    ηv, ρv = value(η_), value(ρ)
    y = Φ(ℓ, ηv, ρv)
    dη = _fd_deriv(x -> Φ(ℓ, x, ρv), ηv)
    dρ = _fd_deriv(x -> Φ(ℓ, ηv, x), ρv)
    return complex(Dual{T}(real(y), real(dη) * partials(η_) + real(dρ) * partials(ρ)), Dual{T}(imag(y), imag(dη) * partials(η_) + imag(dρ) * partials(ρ)))
end

function Φ(ℓ::Dual{T}, η_::Dual{T}, ρ::Dual{T}) where {T}
    ℓv, ηv, ρv = value(ℓ), value(η_), value(ρ)
    y = Φ(ℓv, ηv, ρv)
    dℓ = _fd_deriv(x -> Φ(x, ηv, ρv), ℓv)
    dη = _fd_deriv(x -> Φ(ℓv, x, ρv), ηv)
    dρ = _fd_deriv(x -> Φ(ℓv, ηv, x), ρv)
    return complex(
        Dual{T}(real(y), real(dℓ) * partials(ℓ) + real(dη) * partials(η_) + real(dρ) * partials(ρ)),
        Dual{T}(imag(y), imag(dℓ) * partials(ℓ) + imag(dη) * partials(η_) + imag(dρ) * partials(ρ)),
    )
end

# ── Coulomb: G (real-valued, FD) ───────────────────────────────────────────────

function G(ℓ::Dual{T}, η_::Real, ρ::Real) where {T}
    ℓv = value(ℓ)
    return Dual{T}(G(ℓv, η_, ρ), _fd_deriv(x -> G(x, η_, ρ), ℓv) * partials(ℓ))
end

function G(ℓ::Real, η_::Dual{T}, ρ::Real) where {T}
    ηv = value(η_)
    return Dual{T}(G(ℓ, ηv, ρ), _fd_deriv(x -> G(ℓ, x, ρ), ηv) * partials(η_))
end

function G(ℓ::Real, η_::Real, ρ::Dual{T}) where {T}
    ρv = value(ρ)
    return Dual{T}(G(ℓ, η_, ρv), _fd_deriv(x -> G(ℓ, η_, x), ρv) * partials(ρ))
end

function G(ℓ::Dual{T}, η_::Dual{T}, ρ::Real) where {T}
    ℓv, ηv = value(ℓ), value(η_)
    y = G(ℓv, ηv, ρ)
    return Dual{T}(y, _fd_deriv(x -> G(x, ηv, ρ), ℓv) * partials(ℓ) + _fd_deriv(x -> G(ℓv, x, ρ), ηv) * partials(η_))
end

function G(ℓ::Dual{T}, η_::Real, ρ::Dual{T}) where {T}
    ℓv, ρv = value(ℓ), value(ρ)
    y = G(ℓv, η_, ρv)
    return Dual{T}(y, _fd_deriv(x -> G(x, η_, ρv), ℓv) * partials(ℓ) + _fd_deriv(x -> G(ℓv, η_, x), ρv) * partials(ρ))
end

function G(ℓ::Real, η_::Dual{T}, ρ::Dual{T}) where {T}
    ηv, ρv = value(η_), value(ρ)
    y = G(ℓ, ηv, ρv)
    return Dual{T}(y, _fd_deriv(x -> G(ℓ, x, ρv), ηv) * partials(η_) + _fd_deriv(x -> G(ℓ, ηv, x), ρv) * partials(ρ))
end

function G(ℓ::Dual{T}, η_::Dual{T}, ρ::Dual{T}) where {T}
    ℓv, ηv, ρv = value(ℓ), value(η_), value(ρ)
    y = G(ℓv, ηv, ρv)
    return Dual{T}(
        y,
        _fd_deriv(x -> G(x, ηv, ρv), ℓv) * partials(ℓ) +
            _fd_deriv(x -> G(ℓv, x, ρv), ηv) * partials(η_) +
            _fd_deriv(x -> G(ℓv, ηv, x), ρv) * partials(ρ),
    )
end

# ── Coulomb: M_regularized (real-valued, FD) ───────────────────────────────────

function M_regularized(a::Dual{T}, b::Real, c::Real) where {T}
    av = value(a)
    return Dual{T}(M_regularized(av, b, c), _fd_deriv(x -> M_regularized(x, b, c), av) * partials(a))
end

function M_regularized(a::Real, b::Dual{T}, c::Real) where {T}
    bv = value(b)
    return Dual{T}(M_regularized(a, bv, c), _fd_deriv(x -> M_regularized(a, x, c), bv) * partials(b))
end

function M_regularized(a::Real, b::Real, c::Dual{T}) where {T}
    cv = value(c)
    return Dual{T}(M_regularized(a, b, cv), _fd_deriv(x -> M_regularized(a, b, x), cv) * partials(c))
end

function M_regularized(a::Dual{T}, b::Dual{T}, c::Real) where {T}
    av, bv = value(a), value(b)
    y = M_regularized(av, bv, c)
    return Dual{T}(y, _fd_deriv(x -> M_regularized(x, bv, c), av) * partials(a) + _fd_deriv(x -> M_regularized(av, x, c), bv) * partials(b))
end

function M_regularized(a::Dual{T}, b::Real, c::Dual{T}) where {T}
    av, cv = value(a), value(c)
    y = M_regularized(av, b, cv)
    return Dual{T}(y, _fd_deriv(x -> M_regularized(x, b, cv), av) * partials(a) + _fd_deriv(x -> M_regularized(av, b, x), cv) * partials(c))
end

function M_regularized(a::Real, b::Dual{T}, c::Dual{T}) where {T}
    bv, cv = value(b), value(c)
    y = M_regularized(a, bv, cv)
    return Dual{T}(y, _fd_deriv(x -> M_regularized(a, x, cv), bv) * partials(b) + _fd_deriv(x -> M_regularized(a, bv, x), cv) * partials(c))
end

function M_regularized(a::Dual{T}, b::Dual{T}, c::Dual{T}) where {T}
    av, bv, cv = value(a), value(b), value(c)
    y = M_regularized(av, bv, cv)
    return Dual{T}(
        y,
        _fd_deriv(x -> M_regularized(x, bv, cv), av) * partials(a) +
            _fd_deriv(x -> M_regularized(av, x, cv), bv) * partials(b) +
            _fd_deriv(x -> M_regularized(av, bv, x), cv) * partials(c),
    )
end

# ── Coulomb: w (FD w.r.t. η; ℓ::Integer is not differentiable) ────────────────

function w(ℓ::Integer, η_::Dual{T}) where {T}
    ηv = value(η_)
    return Dual{T}(w(ℓ, ηv), _fd_deriv(x -> w(ℓ, x), ηv) * partials(η_))
end

# ── Fresnel (analytic derivatives) ─────────────────────────────────────────────

function FresnelC(z::Dual{T}) where {T}
    zv = value(z)
    return Dual{T}(FresnelC(zv), cos((π / 2) * zv^2) * partials(z))
end

function FresnelS(z::Dual{T}) where {T}
    zv = value(z)
    return Dual{T}(FresnelS(zv), sin((π / 2) * zv^2) * partials(z))
end

function FresnelE(z::Dual{T}) where {T}
    zv = value(z)
    y = FresnelE(zv)
    dy = exp(im * (π / 2) * zv^2)
    return complex(Dual{T}(real(y), real(dy) * partials(z)), Dual{T}(imag(y), imag(dy) * partials(z)))
end

fresnel(z::Dual{T}) where {T} = (FresnelC(z), FresnelS(z), FresnelE(z))

# ── Dawson integral (analytic derivative) ──────────────────────────────────────

function dawson(x::Dual{T}) where {T}
    xv = value(x)
    y = dawson(xv)
    derivative = muladd(-2 * xv, y, one(xv))
    return Dual{T}(y, derivative * partials(x))
end

# ── Clausen (analytic dCl_n/dθ = ±Cl_{n-1}(θ) for n≥2, FD for n=1) ──────────

function Clausen(n::Int, θv::Dual{T}; N::Int = 10, m::Int = 20) where {T}
    θval = value(θv)
    y = Clausen(n, θval; N = N, m = m)
    dθ = if n == 1
        _fd_deriv(x -> Clausen(1, x; N = N, m = m), θval)
    elseif iseven(n)
        Clausen(n - 1, θval; N = N, m = m)
    else
        -Clausen(n - 1, θval; N = N, m = m)
    end
    return Dual{T}(y, dθ * partials(θv))
end

# ── Bose–Einstein: dBₖ/dη = Bₖ₋₁ ────────────────────────────────────────────────

# Promote before subtraction, including unsigned orders and nested Dual inputs.
function _bose_lower_order(k::Real, x::Real)
    T = float(promote_type(typeof(k), typeof(x)))
    return T(k) - one(T)
end
_bose_lower_order(k::Real, x::Dual) = _bose_lower_order(k, value(x))

# Continue below the lowest tabulated order for derivatives, including nested Duals.
function _bose_continued(k::Real, x::Real)
    k >= -4.5 && return BoseEinsteinIntegralNorm(k, x)
    T = float(promote_type(typeof(k), typeof(x)))
    return FewSpecialFunctions._bose_normalized(T(k), T(x))
end

function _bose_continued(k::Real, x::Dual{T}) where {T}
    xv = value(x)
    return Dual{T}(_bose_continued(k, xv), _bose_continued(_bose_lower_order(k, xv), xv) * partials(x))
end

# Keep Γ(k+1) fixed while lowering the polylogarithm order for derivatives.
function _bose_scaled(k::Real, order::Real, x::Real)
    T = float(promote_type(typeof(k), typeof(order), typeof(x)))
    if T === Float16 || T === Float32
        return T(FewSpecialFunctions._bose_unnormalized(Float64(k), Float64(x), Float64(order)))
    end
    return FewSpecialFunctions._bose_unnormalized(T(k), T(x), T(order))
end

function _bose_scaled(k::Real, order::Real, x::Dual{T}) where {T}
    xv = value(x)
    return Dual{T}(_bose_scaled(k, order, xv), _bose_scaled(k, _bose_lower_order(order, xv), xv) * partials(x))
end

function BoseEinsteinIntegralNorm(k::Real, x::Dual{T}) where {T}
    xv = value(x)
    y = BoseEinsteinIntegralNorm(k, xv)
    return Dual{T}(y, _bose_continued(_bose_lower_order(k, xv), xv) * partials(x))
end

function BoseEinsteinIntegral(k::Real, x::Dual{T}) where {T}
    xv = value(x)
    y = BoseEinsteinIntegral(k, xv)
    return Dual{T}(y, _bose_scaled(k, _bose_lower_order(k, xv), xv) * partials(x))
end

for f in (:BoseEinsteinIntegral, :BoseEinsteinIntegralNorm)
    @eval begin
        $f(k::Dual, x::Real) = throw(DomainError(k, "differentiation in the discrete order k is not supported"))
        $f(k::Dual, x::Dual) = throw(DomainError(k, "differentiation in the discrete order k is not supported"))
    end
end

# ── Fermi-Dirac (FD) ───────────────────────────────────────────────────────────

function FermiDiracIntegral(j::Dual{T}, x::Real) where {T}
    jv = value(j)
    return Dual{T}(FermiDiracIntegral(jv, x), _fd_deriv(t -> FermiDiracIntegral(t, x), jv) * partials(j))
end

function FermiDiracIntegral(j::Real, x::Dual{T}) where {T}
    xv = value(x)
    return Dual{T}(FermiDiracIntegral(j, xv), _fd_deriv(t -> FermiDiracIntegral(j, t), xv) * partials(x))
end

function FermiDiracIntegral(j::Dual{T}, x::Dual{T}) where {T}
    jv, xv = value(j), value(x)
    y = FermiDiracIntegral(jv, xv)
    return Dual{T}(
        y,
        _fd_deriv(t -> FermiDiracIntegral(t, xv), jv) * partials(j) +
            _fd_deriv(t -> FermiDiracIntegral(jv, t), xv) * partials(x),
    )
end

function FermiDiracIntegralNorm(j::Dual{T}, x::Real) where {T}
    jv = value(j)
    return Dual{T}(FermiDiracIntegralNorm(jv, x), _fd_deriv(t -> FermiDiracIntegralNorm(t, x), jv) * partials(j))
end

function FermiDiracIntegralNorm(j::Real, x::Dual{T}) where {T}
    xv = value(x)
    return Dual{T}(FermiDiracIntegralNorm(j, xv), _fd_deriv(t -> FermiDiracIntegralNorm(j, t), xv) * partials(x))
end

function FermiDiracIntegralNorm(j::Dual{T}, x::Dual{T}) where {T}
    jv, xv = value(j), value(x)
    y = FermiDiracIntegralNorm(jv, xv)
    return Dual{T}(
        y,
        _fd_deriv(t -> FermiDiracIntegralNorm(t, xv), jv) * partials(j) +
            _fd_deriv(t -> FermiDiracIntegralNorm(jv, t), xv) * partials(x),
    )
end

# ── MarcumQ (analytic dQ/db = dQdb; FD for M and a) ───────────────────────────

function MarcumQ(M::Dual{T}, a::Real, b::Real) where {T}
    Mv = value(M)
    return Dual{T}(MarcumQ(Mv, a, b), _fd_deriv(x -> MarcumQ(x, a, b), Mv) * partials(M))
end

function MarcumQ(M::Real, a::Dual{T}, b::Real) where {T}
    av = value(a)
    return Dual{T}(MarcumQ(M, av, b), _fd_deriv(x -> MarcumQ(M, x, b), av) * partials(a))
end

function MarcumQ(M::Real, a::Real, b::Dual{T}) where {T}
    bv = value(b)
    return Dual{T}(MarcumQ(M, a, bv), dQdb(M, a, bv) * partials(b))
end

function MarcumQ(M::Dual{T}, a::Dual{T}, b::Real) where {T}
    Mv, av = value(M), value(a)
    y = MarcumQ(Mv, av, b)
    return Dual{T}(
        y,
        _fd_deriv(x -> MarcumQ(x, av, b), Mv) * partials(M) +
            _fd_deriv(x -> MarcumQ(Mv, x, b), av) * partials(a),
    )
end

function MarcumQ(M::Dual{T}, a::Real, b::Dual{T}) where {T}
    Mv, bv = value(M), value(b)
    y = MarcumQ(Mv, a, bv)
    return Dual{T}(y, _fd_deriv(x -> MarcumQ(x, a, bv), Mv) * partials(M) + dQdb(Mv, a, bv) * partials(b))
end

function MarcumQ(M::Real, a::Dual{T}, b::Dual{T}) where {T}
    av, bv = value(a), value(b)
    y = MarcumQ(M, av, bv)
    return Dual{T}(y, _fd_deriv(x -> MarcumQ(M, x, bv), av) * partials(a) + dQdb(M, av, bv) * partials(b))
end

function MarcumQ(M::Dual{T}, a::Dual{T}, b::Dual{T}) where {T}
    Mv, av, bv = value(M), value(a), value(b)
    y = MarcumQ(Mv, av, bv)
    return Dual{T}(
        y,
        _fd_deriv(x -> MarcumQ(x, av, bv), Mv) * partials(M) +
            _fd_deriv(x -> MarcumQ(Mv, x, bv), av) * partials(a) +
            dQdb(Mv, av, bv) * partials(b),
    )
end

# MarcumQ(a, b) — 2-arg convenience form (M=1 fixed)

function MarcumQ(a::Dual{T}, b::Real) where {T}
    av = value(a)
    return Dual{T}(MarcumQ(av, b), _fd_deriv(x -> MarcumQ(x, b), av) * partials(a))
end

function MarcumQ(a::Real, b::Dual{T}) where {T}
    bv = value(b)
    return Dual{T}(MarcumQ(a, bv), dQdb(a, bv) * partials(b))
end

function MarcumQ(a::Dual{T}, b::Dual{T}) where {T}
    av, bv = value(a), value(b)
    y = MarcumQ(av, bv)
    return Dual{T}(y, _fd_deriv(x -> MarcumQ(x, bv), av) * partials(a) + dQdb(av, bv) * partials(b))
end

# ── dQdb (FD for all arguments) ────────────────────────────────────────────────

function dQdb(M::Dual{T}, a::Real, b::Real) where {T}
    Mv = value(M)
    return Dual{T}(dQdb(Mv, a, b), _fd_deriv(x -> dQdb(x, a, b), Mv) * partials(M))
end

function dQdb(M::Real, a::Dual{T}, b::Real) where {T}
    av = value(a)
    return Dual{T}(dQdb(M, av, b), _fd_deriv(x -> dQdb(M, x, b), av) * partials(a))
end

function dQdb(M::Real, a::Real, b::Dual{T}) where {T}
    bv = value(b)
    return Dual{T}(dQdb(M, a, bv), _fd_deriv(x -> dQdb(M, a, x), bv) * partials(b))
end

function dQdb(M::Dual{T}, a::Dual{T}, b::Real) where {T}
    Mv, av = value(M), value(a)
    y = dQdb(Mv, av, b)
    return Dual{T}(y, _fd_deriv(x -> dQdb(x, av, b), Mv) * partials(M) + _fd_deriv(x -> dQdb(Mv, x, b), av) * partials(a))
end

function dQdb(M::Dual{T}, a::Real, b::Dual{T}) where {T}
    Mv, bv = value(M), value(b)
    y = dQdb(Mv, a, bv)
    return Dual{T}(y, _fd_deriv(x -> dQdb(x, a, bv), Mv) * partials(M) + _fd_deriv(x -> dQdb(Mv, a, x), bv) * partials(b))
end

function dQdb(M::Real, a::Dual{T}, b::Dual{T}) where {T}
    av, bv = value(a), value(b)
    y = dQdb(M, av, bv)
    return Dual{T}(y, _fd_deriv(x -> dQdb(M, x, bv), av) * partials(a) + _fd_deriv(x -> dQdb(M, av, x), bv) * partials(b))
end

function dQdb(M::Dual{T}, a::Dual{T}, b::Dual{T}) where {T}
    Mv, av, bv = value(M), value(a), value(b)
    y = dQdb(Mv, av, bv)
    return Dual{T}(
        y,
        _fd_deriv(x -> dQdb(x, av, bv), Mv) * partials(M) +
            _fd_deriv(x -> dQdb(Mv, x, bv), av) * partials(a) +
            _fd_deriv(x -> dQdb(Mv, av, x), bv) * partials(b),
    )
end

# dQdb(a, b) — 2-arg convenience form (M=1 fixed), FD

function dQdb(a::Dual{T}, b::Real) where {T}
    av = value(a)
    return Dual{T}(dQdb(av, b), _fd_deriv(x -> dQdb(x, b), av) * partials(a))
end

function dQdb(a::Real, b::Dual{T}) where {T}
    bv = value(b)
    return Dual{T}(dQdb(a, bv), _fd_deriv(x -> dQdb(a, x), bv) * partials(b))
end

function dQdb(a::Dual{T}, b::Dual{T}) where {T}
    av, bv = value(a), value(b)
    y = dQdb(av, bv)
    return Dual{T}(y, _fd_deriv(x -> dQdb(x, bv), av) * partials(a) + _fd_deriv(x -> dQdb(av, x), bv) * partials(b))
end

# ── Debye function (FD) ────────────────────────────────────────────────────────

# 2-arg form: debye_function(β, x)

function debye_function(β::Dual{T}, x::Real; tol = 1.0e-35, max_terms = 2000) where {T}
    βv = value(β)
    y = debye_function(βv, x; tol = tol, max_terms = max_terms)
    dβ = _fd_deriv(t -> debye_function(t, x; tol = tol, max_terms = max_terms), βv)
    return Dual{T}(y, dβ * partials(β))
end

function debye_function(β::Real, x::Dual{T}; tol = 1.0e-35, max_terms = 2000) where {T}
    xv = value(x)
    y = debye_function(β, xv; tol = tol, max_terms = max_terms)
    dx = _fd_deriv(t -> debye_function(β, t; tol = tol, max_terms = max_terms), xv)
    return Dual{T}(y, dx * partials(x))
end

function debye_function(β::Dual{T}, x::Dual{T}; tol = 1.0e-35, max_terms = 2000) where {T}
    βv, xv = value(β), value(x)
    y = debye_function(βv, xv; tol = tol, max_terms = max_terms)
    dβ = _fd_deriv(t -> debye_function(t, xv; tol = tol, max_terms = max_terms), βv)
    dx = _fd_deriv(t -> debye_function(βv, t; tol = tol, max_terms = max_terms), xv)
    return Dual{T}(y, dβ * partials(β) + dx * partials(x))
end

# 3-arg form: debye_function(n, β, x)

function debye_function(n::Dual{T}, β::Real, x::Real; tol = 1.0e-35, max_terms = 2000) where {T}
    nv = value(n)
    y = debye_function(nv, β, x; tol = tol, max_terms = max_terms)
    dn = _fd_deriv(t -> debye_function(t, β, x; tol = tol, max_terms = max_terms), nv)
    return Dual{T}(y, dn * partials(n))
end

function debye_function(n::Real, β::Dual{T}, x::Real; tol = 1.0e-35, max_terms = 2000) where {T}
    βv = value(β)
    y = debye_function(n, βv, x; tol = tol, max_terms = max_terms)
    dβ = _fd_deriv(t -> debye_function(n, t, x; tol = tol, max_terms = max_terms), βv)
    return Dual{T}(y, dβ * partials(β))
end

function debye_function(n::Real, β::Real, x::Dual{T}; tol = 1.0e-35, max_terms = 2000) where {T}
    xv = value(x)
    y = debye_function(n, β, xv; tol = tol, max_terms = max_terms)
    dx = _fd_deriv(t -> debye_function(n, β, t; tol = tol, max_terms = max_terms), xv)
    return Dual{T}(y, dx * partials(x))
end

function debye_function(n::Dual{T}, β::Dual{T}, x::Real; tol = 1.0e-35, max_terms = 2000) where {T}
    nv, βv = value(n), value(β)
    y = debye_function(nv, βv, x; tol = tol, max_terms = max_terms)
    dn = _fd_deriv(t -> debye_function(t, βv, x; tol = tol, max_terms = max_terms), nv)
    dβ = _fd_deriv(t -> debye_function(nv, t, x; tol = tol, max_terms = max_terms), βv)
    return Dual{T}(y, dn * partials(n) + dβ * partials(β))
end

function debye_function(n::Dual{T}, β::Real, x::Dual{T}; tol = 1.0e-35, max_terms = 2000) where {T}
    nv, xv = value(n), value(x)
    y = debye_function(nv, β, xv; tol = tol, max_terms = max_terms)
    dn = _fd_deriv(t -> debye_function(t, β, xv; tol = tol, max_terms = max_terms), nv)
    dx = _fd_deriv(t -> debye_function(nv, β, t; tol = tol, max_terms = max_terms), xv)
    return Dual{T}(y, dn * partials(n) + dx * partials(x))
end

function debye_function(n::Real, β::Dual{T}, x::Dual{T}; tol = 1.0e-35, max_terms = 2000) where {T}
    βv, xv = value(β), value(x)
    y = debye_function(n, βv, xv; tol = tol, max_terms = max_terms)
    dβ = _fd_deriv(t -> debye_function(n, t, xv; tol = tol, max_terms = max_terms), βv)
    dx = _fd_deriv(t -> debye_function(n, βv, t; tol = tol, max_terms = max_terms), xv)
    return Dual{T}(y, dβ * partials(β) + dx * partials(x))
end

function debye_function(n::Dual{T}, β::Dual{T}, x::Dual{T}; tol = 1.0e-35, max_terms = 2000) where {T}
    nv, βv, xv = value(n), value(β), value(x)
    y = debye_function(nv, βv, xv; tol = tol, max_terms = max_terms)
    dn = _fd_deriv(t -> debye_function(t, βv, xv; tol = tol, max_terms = max_terms), nv)
    dβ = _fd_deriv(t -> debye_function(nv, t, xv; tol = tol, max_terms = max_terms), βv)
    dx = _fd_deriv(t -> debye_function(nv, βv, t; tol = tol, max_terms = max_terms), xv)
    return Dual{T}(y, dn * partials(n) + dβ * partials(β) + dx * partials(x))
end

# Whittaker: analytic argument derivatives, finite differences in κ and μ.
function _whittaker_lift(f, df, args, ::Type{T}) where {T}
    vals = value.(args)
    y = f(vals...)
    p = partials(args[findfirst(x -> x isa Dual, args)])
    pr, pi = zero(p), zero(p)
    for i in 1:3
        args[i] isa Dual || continue
        d = i == 3 ? df(vals...) : _fd_deriv(t -> f(ntuple(j -> j == i ? t : vals[j], 3)...), vals[i])
        pr += real(d) * partials(args[i])
        pi += imag(d) * partials(args[i])
    end
    return y isa Real ? Dual{T}(y, pr) : complex(Dual{T}(real(y), pr), Dual{T}(imag(y), pi))
end

_whittaker_second(f, κ, μ, z) = (1 / 4 - κ / z - (1 / 4 - μ^2) / z^2) * f(κ, μ, z)

for (f, df) in (
        (:WhittakerM, :dWhittakerM), (:WhittakerW, :dWhittakerW),
        (:dWhittakerM, :((κ, μ, z) -> _whittaker_second(WhittakerM, κ, μ, z))),
        (:dWhittakerW, :((κ, μ, z) -> _whittaker_second(WhittakerW, κ, μ, z))),
    )
    for mask in 1:7
        args = [:($(name)::$(iszero(mask & (1 << (i - 1))) ? :Number : :(Dual{T}))) for (i, name) in enumerate((:κ, :μ, :z))]
        @eval $f($(args...)) where {T} = _whittaker_lift($f, $df, (κ, μ, z), T)
    end
end

function _cylinder_scaled_dx(a, x, kind)
    af, xf = promote(float(a), float(x))
    y, dy = FewSpecialFunctions._cylinder_scaled_pair(af, xf, kind)
    dlogscale = sqrt(max(zero(xf), xf^2 / 4 + af))
    return dy + (kind === :U ? dlogscale : -dlogscale) * y
end

for (f, kind) in ((:U_scaled, :U), (:V_scaled, :V))
    @eval begin
        function $f(a::Dual{T}, x::Real) where {T}
            av = value(a)
            return Dual{T}($f(av, x), _fd_deriv(t -> $f(t, x), av) * partials(a))
        end
        function $f(a::Real, x::Dual{T}) where {T}
            xv = value(x)
            return Dual{T}($f(a, xv), _cylinder_scaled_dx(a, xv, $(QuoteNode(kind))) * partials(x))
        end
        function $f(a::Dual{T}, x::Dual{T}) where {T}
            av, xv = value(a), value(x)
            return Dual{T}($f(av, xv), _fd_deriv(t -> $f(t, xv), av) * partials(a) + _cylinder_scaled_dx(av, xv, $(QuoteNode(kind))) * partials(x))
        end
    end
end

# ── Parabolic cylinder: U, V, W ────────────────────────────────────────────────
# ∂U/∂x = dU(a, x),  ∂V/∂x = dV(a, x),  ∂W/∂x = dW(a, x)

function U(a::Dual{T}, x::Real) where {T}
    av = value(a)
    return Dual{T}(U(av, x), _fd_deriv(t -> U(t, x), av) * partials(a))
end

function U(a::Real, x::Dual{T}) where {T}
    xv = value(x)
    return Dual{T}(U(a, xv), dU(a, xv) * partials(x))
end

function U(a::Dual{T}, x::Dual{T}) where {T}
    av, xv = value(a), value(x)
    return Dual{T}(U(av, xv), _fd_deriv(t -> U(t, xv), av) * partials(a) + dU(av, xv) * partials(x))
end

function V(a::Dual{T}, x::Real) where {T}
    av = value(a)
    return Dual{T}(V(av, x), _fd_deriv(t -> V(t, x), av) * partials(a))
end

function V(a::Real, x::Dual{T}) where {T}
    xv = value(x)
    return Dual{T}(V(a, xv), dV(a, xv) * partials(x))
end

function V(a::Dual{T}, x::Dual{T}) where {T}
    av, xv = value(a), value(x)
    return Dual{T}(V(av, xv), _fd_deriv(t -> V(t, xv), av) * partials(a) + dV(av, xv) * partials(x))
end

function W(a::Dual{T}, x::Real) where {T}
    av = value(a)
    return Dual{T}(W(av, x), _fd_deriv(t -> W(t, x), av) * partials(a))
end

function W(a::Real, x::Dual{T}) where {T}
    xv = value(x)
    return Dual{T}(W(a, xv), dW(a, xv) * partials(x))
end

function W(a::Dual{T}, x::Dual{T}) where {T}
    av, xv = value(a), value(x)
    return Dual{T}(W(av, xv), _fd_deriv(t -> W(t, xv), av) * partials(a) + dW(av, xv) * partials(x))
end

# ── Parabolic cylinder: dU, dV ─────────────────────────────────────────────────
# d(dU)/dx = (x²/4 + a) * U(a, x)  from the parabolic cylinder ODE U'' = (x²/4 + a)U

function dU(a::Dual{T}, x::Real) where {T}
    av = value(a)
    return Dual{T}(dU(av, x), _fd_deriv(t -> dU(t, x), av) * partials(a))
end

function dU(a::Real, x::Dual{T}) where {T}
    xv = value(x)
    return Dual{T}(dU(a, xv), (xv^2 / 4 + a) * U(a, xv) * partials(x))
end

function dU(a::Dual{T}, x::Dual{T}) where {T}
    av, xv = value(a), value(x)
    return Dual{T}(dU(av, xv), _fd_deriv(t -> dU(t, xv), av) * partials(a) + (xv^2 / 4 + av) * U(av, xv) * partials(x))
end

function dV(a::Dual{T}, x::Real) where {T}
    av = value(a)
    return Dual{T}(dV(av, x), _fd_deriv(t -> dV(t, x), av) * partials(a))
end

function dV(a::Real, x::Dual{T}) where {T}
    xv = value(x)
    return Dual{T}(dV(a, xv), (xv^2 / 4 + a) * V(a, xv) * partials(x))
end

function dV(a::Dual{T}, x::Dual{T}) where {T}
    av, xv = value(a), value(x)
    return Dual{T}(dV(av, xv), _fd_deriv(t -> dV(t, xv), av) * partials(a) + (xv^2 / 4 + av) * V(av, xv) * partials(x))
end

# ── Parabolic cylinder: dW (W'' = (a - x²/4)W) ────────────────────────────────

function dW(a::Dual{T}, x::Real) where {T}
    av = value(a)
    return Dual{T}(dW(av, x), _fd_deriv(t -> dW(t, x), av) * partials(a))
end

function dW(a::Real, x::Dual{T}) where {T}
    xv = value(x)
    return Dual{T}(dW(a, xv), (a - xv^2 / 4) * W(a, xv) * partials(x))
end

function dW(a::Dual{T}, x::Dual{T}) where {T}
    av, xv = value(a), value(x)
    return Dual{T}(dW(av, xv), _fd_deriv(t -> dW(t, xv), av) * partials(a) + (av - xv^2 / 4) * W(av, xv) * partials(x))
end

end
