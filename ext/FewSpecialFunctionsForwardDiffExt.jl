module FewSpecialFunctionsForwardDiffExt

using FewSpecialFunctions
using ForwardDiff
using ForwardDiff: Dual, partials, value

import FewSpecialFunctions:
    η, C, θ, F, D⁺, D⁻, H⁺, H⁻, F_imag, G, M_regularized, Φ, w,
    debye_function,
    fresnel, FresnelS, FresnelC, FresnelE,
    Clausen,
    FermiDiracIntegral, FermiDiracIntegralNorm,
    MarcumQ, dQdb,
    U, V, W, dU, dV, dW

_fd_step(x::T) where {T <: Real} = cbrt(eps(T)) * (abs(x) + one(T))
_fd_deriv(f, x) = (h = _fd_step(x); (f(x + h) - f(x - h)) / (2h))

# For complex-valued functions: split the derivative into real and imaginary Dual parts.
_complex_dual(::Type{T}, y, dy, p) where {T} =
    complex(Dual{T}(real(y), real(dy) * p), Dual{T}(imag(y), imag(dy) * p))

# ── Coulomb: η ─────────────────────────────────────────────────────────────────
# dη/dϵ = -1 / (2 ϵ^(3/2))

function η(ϵ::Dual{T}) where {T}
    ϵv = value(ϵ)
    Dual{T}(η(ϵv), -one(ϵv) / (2 * ϵv * sqrt(ϵv)) * partials(ϵ))
end

# dη/da = -1/(a²k),  dη/dk = -1/(ak²)

function η(a::Dual{T}, k::Real) where {T}
    av = value(a)
    Dual{T}(η(av, k), -one(av) / (av^2 * k) * partials(a))
end

function η(a::Real, k::Dual{T}) where {T}
    kv = value(k)
    Dual{T}(η(a, kv), -one(kv) / (a * kv^2) * partials(k))
end

function η(a::Dual{T}, k::Dual{T}) where {T}
    av, kv = value(a), value(k)
    Dual{T}(η(av, kv), (-one(av) / (av^2 * kv)) * partials(a) + (-one(kv) / (av * kv^2)) * partials(k))
end

# ── Coulomb: C (real-valued, FD) ───────────────────────────────────────────────

function C(ℓ::Dual{T}, η_::Real) where {T}
    ℓv = value(ℓ)
    Dual{T}(C(ℓv, η_), _fd_deriv(x -> C(x, η_), ℓv) * partials(ℓ))
end

function C(ℓ::Real, η_::Dual{T}) where {T}
    ηv = value(η_)
    Dual{T}(C(ℓ, ηv), _fd_deriv(x -> C(ℓ, x), ηv) * partials(η_))
end

function C(ℓ::Dual{T}, η_::Dual{T}) where {T}
    ℓv, ηv = value(ℓ), value(η_)
    y = C(ℓv, ηv)
    Dual{T}(y, _fd_deriv(x -> C(x, ηv), ℓv) * partials(ℓ) + _fd_deriv(x -> C(ℓv, x), ηv) * partials(η_))
end

# ── Coulomb: D⁺, D⁻ (complex-valued, FD) ──────────────────────────────────────

function D⁺(ℓ::Dual{T}, η_::Real) where {T}
    ℓv = value(ℓ)
    y = D⁺(ℓv, η_)
    _complex_dual(T, y, _fd_deriv(x -> D⁺(x, η_), ℓv), partials(ℓ))
end

function D⁺(ℓ::Real, η_::Dual{T}) where {T}
    ηv = value(η_)
    y = D⁺(ℓ, ηv)
    _complex_dual(T, y, _fd_deriv(x -> D⁺(ℓ, x), ηv), partials(η_))
end

function D⁺(ℓ::Dual{T}, η_::Dual{T}) where {T}
    ℓv, ηv = value(ℓ), value(η_)
    y = D⁺(ℓv, ηv)
    dℓ = _fd_deriv(x -> D⁺(x, ηv), ℓv)
    dη = _fd_deriv(x -> D⁺(ℓv, x), ηv)
    complex(
        Dual{T}(real(y), real(dℓ) * partials(ℓ) + real(dη) * partials(η_)),
        Dual{T}(imag(y), imag(dℓ) * partials(ℓ) + imag(dη) * partials(η_)),
    )
end

function D⁻(ℓ::Dual{T}, η_::Real) where {T}
    ℓv = value(ℓ)
    y = D⁻(ℓv, η_)
    _complex_dual(T, y, _fd_deriv(x -> D⁻(x, η_), ℓv), partials(ℓ))
end

function D⁻(ℓ::Real, η_::Dual{T}) where {T}
    ηv = value(η_)
    y = D⁻(ℓ, ηv)
    _complex_dual(T, y, _fd_deriv(x -> D⁻(ℓ, x), ηv), partials(η_))
end

function D⁻(ℓ::Dual{T}, η_::Dual{T}) where {T}
    ℓv, ηv = value(ℓ), value(η_)
    y = D⁻(ℓv, ηv)
    dℓ = _fd_deriv(x -> D⁻(x, ηv), ℓv)
    dη = _fd_deriv(x -> D⁻(ℓv, x), ηv)
    complex(
        Dual{T}(real(y), real(dℓ) * partials(ℓ) + real(dη) * partials(η_)),
        Dual{T}(imag(y), imag(dℓ) * partials(ℓ) + imag(dη) * partials(η_)),
    )
end

# ── Coulomb: θ (real-valued; analytic dθ/dρ = 1 - η/ρ, FD for ℓ and η) ───────

function θ(ℓ::Dual{T}, η_::Real, ρ::Real) where {T}
    ℓv = value(ℓ)
    Dual{T}(θ(ℓv, η_, ρ), _fd_deriv(x -> θ(x, η_, ρ), ℓv) * partials(ℓ))
end

function θ(ℓ::Real, η_::Dual{T}, ρ::Real) where {T}
    ηv = value(η_)
    Dual{T}(θ(ℓ, ηv, ρ), _fd_deriv(x -> θ(ℓ, x, ρ), ηv) * partials(η_))
end

function θ(ℓ::Real, η_::Real, ρ::Dual{T}) where {T}
    ρv = value(ρ)
    Dual{T}(θ(ℓ, η_, ρv), (one(ρv) - η_ / ρv) * partials(ρ))
end

function θ(ℓ::Dual{T}, η_::Dual{T}, ρ::Real) where {T}
    ℓv, ηv = value(ℓ), value(η_)
    y = θ(ℓv, ηv, ρ)
    Dual{T}(y, _fd_deriv(x -> θ(x, ηv, ρ), ℓv) * partials(ℓ) + _fd_deriv(x -> θ(ℓv, x, ρ), ηv) * partials(η_))
end

function θ(ℓ::Dual{T}, η_::Real, ρ::Dual{T}) where {T}
    ℓv, ρv = value(ℓ), value(ρ)
    y = θ(ℓv, η_, ρv)
    Dual{T}(y, _fd_deriv(x -> θ(x, η_, ρv), ℓv) * partials(ℓ) + (one(ρv) - η_ / ρv) * partials(ρ))
end

function θ(ℓ::Real, η_::Dual{T}, ρ::Dual{T}) where {T}
    ηv, ρv = value(η_), value(ρ)
    y = θ(ℓ, ηv, ρv)
    Dual{T}(y, _fd_deriv(x -> θ(ℓ, x, ρv), ηv) * partials(η_) + (one(ρv) - ηv / ρv) * partials(ρ))
end

function θ(ℓ::Dual{T}, η_::Dual{T}, ρ::Dual{T}) where {T}
    ℓv, ηv, ρv = value(ℓ), value(η_), value(ρ)
    y = θ(ℓv, ηv, ρv)
    Dual{T}(
        y,
        _fd_deriv(x -> θ(x, ηv, ρv), ℓv) * partials(ℓ) +
            _fd_deriv(x -> θ(ℓv, x, ρv), ηv) * partials(η_) +
            (one(ρv) - ηv / ρv) * partials(ρ),
    )
end

# ── Coulomb: F (real-valued, FD) ───────────────────────────────────────────────

function F(ℓ::Dual{T}, η_::Real, ρ::Real) where {T}
    ℓv = value(ℓ)
    Dual{T}(F(ℓv, η_, ρ), _fd_deriv(x -> F(x, η_, ρ), ℓv) * partials(ℓ))
end

function F(ℓ::Real, η_::Dual{T}, ρ::Real) where {T}
    ηv = value(η_)
    Dual{T}(F(ℓ, ηv, ρ), _fd_deriv(x -> F(ℓ, x, ρ), ηv) * partials(η_))
end

function F(ℓ::Real, η_::Real, ρ::Dual{T}) where {T}
    ρv = value(ρ)
    Dual{T}(F(ℓ, η_, ρv), _fd_deriv(x -> F(ℓ, η_, x), ρv) * partials(ρ))
end

function F(ℓ::Dual{T}, η_::Dual{T}, ρ::Real) where {T}
    ℓv, ηv = value(ℓ), value(η_)
    y = F(ℓv, ηv, ρ)
    Dual{T}(y, _fd_deriv(x -> F(x, ηv, ρ), ℓv) * partials(ℓ) + _fd_deriv(x -> F(ℓv, x, ρ), ηv) * partials(η_))
end

function F(ℓ::Dual{T}, η_::Real, ρ::Dual{T}) where {T}
    ℓv, ρv = value(ℓ), value(ρ)
    y = F(ℓv, η_, ρv)
    Dual{T}(y, _fd_deriv(x -> F(x, η_, ρv), ℓv) * partials(ℓ) + _fd_deriv(x -> F(ℓv, η_, x), ρv) * partials(ρ))
end

function F(ℓ::Real, η_::Dual{T}, ρ::Dual{T}) where {T}
    ηv, ρv = value(η_), value(ρ)
    y = F(ℓ, ηv, ρv)
    Dual{T}(y, _fd_deriv(x -> F(ℓ, x, ρv), ηv) * partials(η_) + _fd_deriv(x -> F(ℓ, ηv, x), ρv) * partials(ρ))
end

function F(ℓ::Dual{T}, η_::Dual{T}, ρ::Dual{T}) where {T}
    ℓv, ηv, ρv = value(ℓ), value(η_), value(ρ)
    y = F(ℓv, ηv, ρv)
    Dual{T}(
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
    _complex_dual(T, y, _fd_deriv(x -> H⁺(x, η_, ρ), ℓv), partials(ℓ))
end

function H⁺(ℓ::Real, η_::Dual{T}, ρ::Real) where {T}
    ηv = value(η_)
    y = H⁺(ℓ, ηv, ρ)
    _complex_dual(T, y, _fd_deriv(x -> H⁺(ℓ, x, ρ), ηv), partials(η_))
end

function H⁺(ℓ::Real, η_::Real, ρ::Dual{T}) where {T}
    ρv = value(ρ)
    y = H⁺(ℓ, η_, ρv)
    _complex_dual(T, y, _fd_deriv(x -> H⁺(ℓ, η_, x), ρv), partials(ρ))
end

function H⁺(ℓ::Dual{T}, η_::Dual{T}, ρ::Real) where {T}
    ℓv, ηv = value(ℓ), value(η_)
    y = H⁺(ℓv, ηv, ρ)
    dℓ = _fd_deriv(x -> H⁺(x, ηv, ρ), ℓv)
    dη = _fd_deriv(x -> H⁺(ℓv, x, ρ), ηv)
    complex(Dual{T}(real(y), real(dℓ) * partials(ℓ) + real(dη) * partials(η_)), Dual{T}(imag(y), imag(dℓ) * partials(ℓ) + imag(dη) * partials(η_)))
end

function H⁺(ℓ::Dual{T}, η_::Real, ρ::Dual{T}) where {T}
    ℓv, ρv = value(ℓ), value(ρ)
    y = H⁺(ℓv, η_, ρv)
    dℓ = _fd_deriv(x -> H⁺(x, η_, ρv), ℓv)
    dρ = _fd_deriv(x -> H⁺(ℓv, η_, x), ρv)
    complex(Dual{T}(real(y), real(dℓ) * partials(ℓ) + real(dρ) * partials(ρ)), Dual{T}(imag(y), imag(dℓ) * partials(ℓ) + imag(dρ) * partials(ρ)))
end

function H⁺(ℓ::Real, η_::Dual{T}, ρ::Dual{T}) where {T}
    ηv, ρv = value(η_), value(ρ)
    y = H⁺(ℓ, ηv, ρv)
    dη = _fd_deriv(x -> H⁺(ℓ, x, ρv), ηv)
    dρ = _fd_deriv(x -> H⁺(ℓ, ηv, x), ρv)
    complex(Dual{T}(real(y), real(dη) * partials(η_) + real(dρ) * partials(ρ)), Dual{T}(imag(y), imag(dη) * partials(η_) + imag(dρ) * partials(ρ)))
end

function H⁺(ℓ::Dual{T}, η_::Dual{T}, ρ::Dual{T}) where {T}
    ℓv, ηv, ρv = value(ℓ), value(η_), value(ρ)
    y = H⁺(ℓv, ηv, ρv)
    dℓ = _fd_deriv(x -> H⁺(x, ηv, ρv), ℓv)
    dη = _fd_deriv(x -> H⁺(ℓv, x, ρv), ηv)
    dρ = _fd_deriv(x -> H⁺(ℓv, ηv, x), ρv)
    complex(
        Dual{T}(real(y), real(dℓ) * partials(ℓ) + real(dη) * partials(η_) + real(dρ) * partials(ρ)),
        Dual{T}(imag(y), imag(dℓ) * partials(ℓ) + imag(dη) * partials(η_) + imag(dρ) * partials(ρ)),
    )
end

# ── Coulomb: H⁻ (complex-valued, FD) ──────────────────────────────────────────

function H⁻(ℓ::Dual{T}, η_::Real, ρ::Real) where {T}
    ℓv = value(ℓ)
    y = H⁻(ℓv, η_, ρ)
    _complex_dual(T, y, _fd_deriv(x -> H⁻(x, η_, ρ), ℓv), partials(ℓ))
end

function H⁻(ℓ::Real, η_::Dual{T}, ρ::Real) where {T}
    ηv = value(η_)
    y = H⁻(ℓ, ηv, ρ)
    _complex_dual(T, y, _fd_deriv(x -> H⁻(ℓ, x, ρ), ηv), partials(η_))
end

function H⁻(ℓ::Real, η_::Real, ρ::Dual{T}) where {T}
    ρv = value(ρ)
    y = H⁻(ℓ, η_, ρv)
    _complex_dual(T, y, _fd_deriv(x -> H⁻(ℓ, η_, x), ρv), partials(ρ))
end

function H⁻(ℓ::Dual{T}, η_::Dual{T}, ρ::Real) where {T}
    ℓv, ηv = value(ℓ), value(η_)
    y = H⁻(ℓv, ηv, ρ)
    dℓ = _fd_deriv(x -> H⁻(x, ηv, ρ), ℓv)
    dη = _fd_deriv(x -> H⁻(ℓv, x, ρ), ηv)
    complex(Dual{T}(real(y), real(dℓ) * partials(ℓ) + real(dη) * partials(η_)), Dual{T}(imag(y), imag(dℓ) * partials(ℓ) + imag(dη) * partials(η_)))
end

function H⁻(ℓ::Dual{T}, η_::Real, ρ::Dual{T}) where {T}
    ℓv, ρv = value(ℓ), value(ρ)
    y = H⁻(ℓv, η_, ρv)
    dℓ = _fd_deriv(x -> H⁻(x, η_, ρv), ℓv)
    dρ = _fd_deriv(x -> H⁻(ℓv, η_, x), ρv)
    complex(Dual{T}(real(y), real(dℓ) * partials(ℓ) + real(dρ) * partials(ρ)), Dual{T}(imag(y), imag(dℓ) * partials(ℓ) + imag(dρ) * partials(ρ)))
end

function H⁻(ℓ::Real, η_::Dual{T}, ρ::Dual{T}) where {T}
    ηv, ρv = value(η_), value(ρ)
    y = H⁻(ℓ, ηv, ρv)
    dη = _fd_deriv(x -> H⁻(ℓ, x, ρv), ηv)
    dρ = _fd_deriv(x -> H⁻(ℓ, ηv, x), ρv)
    complex(Dual{T}(real(y), real(dη) * partials(η_) + real(dρ) * partials(ρ)), Dual{T}(imag(y), imag(dη) * partials(η_) + imag(dρ) * partials(ρ)))
end

function H⁻(ℓ::Dual{T}, η_::Dual{T}, ρ::Dual{T}) where {T}
    ℓv, ηv, ρv = value(ℓ), value(η_), value(ρ)
    y = H⁻(ℓv, ηv, ρv)
    dℓ = _fd_deriv(x -> H⁻(x, ηv, ρv), ℓv)
    dη = _fd_deriv(x -> H⁻(ℓv, x, ρv), ηv)
    dρ = _fd_deriv(x -> H⁻(ℓv, ηv, x), ρv)
    complex(
        Dual{T}(real(y), real(dℓ) * partials(ℓ) + real(dη) * partials(η_) + real(dρ) * partials(ρ)),
        Dual{T}(imag(y), imag(dℓ) * partials(ℓ) + imag(dη) * partials(η_) + imag(dρ) * partials(ρ)),
    )
end

# ── Coulomb: F_imag (complex-valued, FD) ───────────────────────────────────────

function F_imag(ℓ::Dual{T}, η_::Real, ρ::Real) where {T}
    ℓv = value(ℓ)
    y = F_imag(ℓv, η_, ρ)
    _complex_dual(T, y, _fd_deriv(x -> F_imag(x, η_, ρ), ℓv), partials(ℓ))
end

function F_imag(ℓ::Real, η_::Dual{T}, ρ::Real) where {T}
    ηv = value(η_)
    y = F_imag(ℓ, ηv, ρ)
    _complex_dual(T, y, _fd_deriv(x -> F_imag(ℓ, x, ρ), ηv), partials(η_))
end

function F_imag(ℓ::Real, η_::Real, ρ::Dual{T}) where {T}
    ρv = value(ρ)
    y = F_imag(ℓ, η_, ρv)
    _complex_dual(T, y, _fd_deriv(x -> F_imag(ℓ, η_, x), ρv), partials(ρ))
end

function F_imag(ℓ::Dual{T}, η_::Dual{T}, ρ::Real) where {T}
    ℓv, ηv = value(ℓ), value(η_)
    y = F_imag(ℓv, ηv, ρ)
    dℓ = _fd_deriv(x -> F_imag(x, ηv, ρ), ℓv)
    dη = _fd_deriv(x -> F_imag(ℓv, x, ρ), ηv)
    complex(Dual{T}(real(y), real(dℓ) * partials(ℓ) + real(dη) * partials(η_)), Dual{T}(imag(y), imag(dℓ) * partials(ℓ) + imag(dη) * partials(η_)))
end

function F_imag(ℓ::Dual{T}, η_::Real, ρ::Dual{T}) where {T}
    ℓv, ρv = value(ℓ), value(ρ)
    y = F_imag(ℓv, η_, ρv)
    dℓ = _fd_deriv(x -> F_imag(x, η_, ρv), ℓv)
    dρ = _fd_deriv(x -> F_imag(ℓv, η_, x), ρv)
    complex(Dual{T}(real(y), real(dℓ) * partials(ℓ) + real(dρ) * partials(ρ)), Dual{T}(imag(y), imag(dℓ) * partials(ℓ) + imag(dρ) * partials(ρ)))
end

function F_imag(ℓ::Real, η_::Dual{T}, ρ::Dual{T}) where {T}
    ηv, ρv = value(η_), value(ρ)
    y = F_imag(ℓ, ηv, ρv)
    dη = _fd_deriv(x -> F_imag(ℓ, x, ρv), ηv)
    dρ = _fd_deriv(x -> F_imag(ℓ, ηv, x), ρv)
    complex(Dual{T}(real(y), real(dη) * partials(η_) + real(dρ) * partials(ρ)), Dual{T}(imag(y), imag(dη) * partials(η_) + imag(dρ) * partials(ρ)))
end

function F_imag(ℓ::Dual{T}, η_::Dual{T}, ρ::Dual{T}) where {T}
    ℓv, ηv, ρv = value(ℓ), value(η_), value(ρ)
    y = F_imag(ℓv, ηv, ρv)
    dℓ = _fd_deriv(x -> F_imag(x, ηv, ρv), ℓv)
    dη = _fd_deriv(x -> F_imag(ℓv, x, ρv), ηv)
    dρ = _fd_deriv(x -> F_imag(ℓv, ηv, x), ρv)
    complex(
        Dual{T}(real(y), real(dℓ) * partials(ℓ) + real(dη) * partials(η_) + real(dρ) * partials(ρ)),
        Dual{T}(imag(y), imag(dℓ) * partials(ℓ) + imag(dη) * partials(η_) + imag(dρ) * partials(ρ)),
    )
end

# ── Coulomb: Φ (complex-valued, FD) ────────────────────────────────────────────

function Φ(ℓ::Dual{T}, η_::Real, ρ::Real) where {T}
    ℓv = value(ℓ)
    y = Φ(ℓv, η_, ρ)
    _complex_dual(T, y, _fd_deriv(x -> Φ(x, η_, ρ), ℓv), partials(ℓ))
end

function Φ(ℓ::Real, η_::Dual{T}, ρ::Real) where {T}
    ηv = value(η_)
    y = Φ(ℓ, ηv, ρ)
    _complex_dual(T, y, _fd_deriv(x -> Φ(ℓ, x, ρ), ηv), partials(η_))
end

function Φ(ℓ::Real, η_::Real, ρ::Dual{T}) where {T}
    ρv = value(ρ)
    y = Φ(ℓ, η_, ρv)
    _complex_dual(T, y, _fd_deriv(x -> Φ(ℓ, η_, x), ρv), partials(ρ))
end

function Φ(ℓ::Dual{T}, η_::Dual{T}, ρ::Real) where {T}
    ℓv, ηv = value(ℓ), value(η_)
    y = Φ(ℓv, ηv, ρ)
    dℓ = _fd_deriv(x -> Φ(x, ηv, ρ), ℓv)
    dη = _fd_deriv(x -> Φ(ℓv, x, ρ), ηv)
    complex(Dual{T}(real(y), real(dℓ) * partials(ℓ) + real(dη) * partials(η_)), Dual{T}(imag(y), imag(dℓ) * partials(ℓ) + imag(dη) * partials(η_)))
end

function Φ(ℓ::Dual{T}, η_::Real, ρ::Dual{T}) where {T}
    ℓv, ρv = value(ℓ), value(ρ)
    y = Φ(ℓv, η_, ρv)
    dℓ = _fd_deriv(x -> Φ(x, η_, ρv), ℓv)
    dρ = _fd_deriv(x -> Φ(ℓv, η_, x), ρv)
    complex(Dual{T}(real(y), real(dℓ) * partials(ℓ) + real(dρ) * partials(ρ)), Dual{T}(imag(y), imag(dℓ) * partials(ℓ) + imag(dρ) * partials(ρ)))
end

function Φ(ℓ::Real, η_::Dual{T}, ρ::Dual{T}) where {T}
    ηv, ρv = value(η_), value(ρ)
    y = Φ(ℓ, ηv, ρv)
    dη = _fd_deriv(x -> Φ(ℓ, x, ρv), ηv)
    dρ = _fd_deriv(x -> Φ(ℓ, ηv, x), ρv)
    complex(Dual{T}(real(y), real(dη) * partials(η_) + real(dρ) * partials(ρ)), Dual{T}(imag(y), imag(dη) * partials(η_) + imag(dρ) * partials(ρ)))
end

function Φ(ℓ::Dual{T}, η_::Dual{T}, ρ::Dual{T}) where {T}
    ℓv, ηv, ρv = value(ℓ), value(η_), value(ρ)
    y = Φ(ℓv, ηv, ρv)
    dℓ = _fd_deriv(x -> Φ(x, ηv, ρv), ℓv)
    dη = _fd_deriv(x -> Φ(ℓv, x, ρv), ηv)
    dρ = _fd_deriv(x -> Φ(ℓv, ηv, x), ρv)
    complex(
        Dual{T}(real(y), real(dℓ) * partials(ℓ) + real(dη) * partials(η_) + real(dρ) * partials(ρ)),
        Dual{T}(imag(y), imag(dℓ) * partials(ℓ) + imag(dη) * partials(η_) + imag(dρ) * partials(ρ)),
    )
end

# ── Coulomb: G (real-valued, FD) ───────────────────────────────────────────────

function G(ℓ::Dual{T}, η_::Real, ρ::Real) where {T}
    ℓv = value(ℓ)
    Dual{T}(G(ℓv, η_, ρ), _fd_deriv(x -> G(x, η_, ρ), ℓv) * partials(ℓ))
end

function G(ℓ::Real, η_::Dual{T}, ρ::Real) where {T}
    ηv = value(η_)
    Dual{T}(G(ℓ, ηv, ρ), _fd_deriv(x -> G(ℓ, x, ρ), ηv) * partials(η_))
end

function G(ℓ::Real, η_::Real, ρ::Dual{T}) where {T}
    ρv = value(ρ)
    Dual{T}(G(ℓ, η_, ρv), _fd_deriv(x -> G(ℓ, η_, x), ρv) * partials(ρ))
end

function G(ℓ::Dual{T}, η_::Dual{T}, ρ::Real) where {T}
    ℓv, ηv = value(ℓ), value(η_)
    y = G(ℓv, ηv, ρ)
    Dual{T}(y, _fd_deriv(x -> G(x, ηv, ρ), ℓv) * partials(ℓ) + _fd_deriv(x -> G(ℓv, x, ρ), ηv) * partials(η_))
end

function G(ℓ::Dual{T}, η_::Real, ρ::Dual{T}) where {T}
    ℓv, ρv = value(ℓ), value(ρ)
    y = G(ℓv, η_, ρv)
    Dual{T}(y, _fd_deriv(x -> G(x, η_, ρv), ℓv) * partials(ℓ) + _fd_deriv(x -> G(ℓv, η_, x), ρv) * partials(ρ))
end

function G(ℓ::Real, η_::Dual{T}, ρ::Dual{T}) where {T}
    ηv, ρv = value(η_), value(ρ)
    y = G(ℓ, ηv, ρv)
    Dual{T}(y, _fd_deriv(x -> G(ℓ, x, ρv), ηv) * partials(η_) + _fd_deriv(x -> G(ℓ, ηv, x), ρv) * partials(ρ))
end

function G(ℓ::Dual{T}, η_::Dual{T}, ρ::Dual{T}) where {T}
    ℓv, ηv, ρv = value(ℓ), value(η_), value(ρ)
    y = G(ℓv, ηv, ρv)
    Dual{T}(
        y,
        _fd_deriv(x -> G(x, ηv, ρv), ℓv) * partials(ℓ) +
            _fd_deriv(x -> G(ℓv, x, ρv), ηv) * partials(η_) +
            _fd_deriv(x -> G(ℓv, ηv, x), ρv) * partials(ρ),
    )
end

# ── Coulomb: M_regularized (real-valued, FD) ───────────────────────────────────

function M_regularized(a::Dual{T}, b::Real, c::Real) where {T}
    av = value(a)
    Dual{T}(M_regularized(av, b, c), _fd_deriv(x -> M_regularized(x, b, c), av) * partials(a))
end

function M_regularized(a::Real, b::Dual{T}, c::Real) where {T}
    bv = value(b)
    Dual{T}(M_regularized(a, bv, c), _fd_deriv(x -> M_regularized(a, x, c), bv) * partials(b))
end

function M_regularized(a::Real, b::Real, c::Dual{T}) where {T}
    cv = value(c)
    Dual{T}(M_regularized(a, b, cv), _fd_deriv(x -> M_regularized(a, b, x), cv) * partials(c))
end

function M_regularized(a::Dual{T}, b::Dual{T}, c::Real) where {T}
    av, bv = value(a), value(b)
    y = M_regularized(av, bv, c)
    Dual{T}(y, _fd_deriv(x -> M_regularized(x, bv, c), av) * partials(a) + _fd_deriv(x -> M_regularized(av, x, c), bv) * partials(b))
end

function M_regularized(a::Dual{T}, b::Real, c::Dual{T}) where {T}
    av, cv = value(a), value(c)
    y = M_regularized(av, b, cv)
    Dual{T}(y, _fd_deriv(x -> M_regularized(x, b, cv), av) * partials(a) + _fd_deriv(x -> M_regularized(av, b, x), cv) * partials(c))
end

function M_regularized(a::Real, b::Dual{T}, c::Dual{T}) where {T}
    bv, cv = value(b), value(c)
    y = M_regularized(a, bv, cv)
    Dual{T}(y, _fd_deriv(x -> M_regularized(a, x, cv), bv) * partials(b) + _fd_deriv(x -> M_regularized(a, bv, x), cv) * partials(c))
end

function M_regularized(a::Dual{T}, b::Dual{T}, c::Dual{T}) where {T}
    av, bv, cv = value(a), value(b), value(c)
    y = M_regularized(av, bv, cv)
    Dual{T}(
        y,
        _fd_deriv(x -> M_regularized(x, bv, cv), av) * partials(a) +
            _fd_deriv(x -> M_regularized(av, x, cv), bv) * partials(b) +
            _fd_deriv(x -> M_regularized(av, bv, x), cv) * partials(c),
    )
end

# ── Coulomb: w (FD w.r.t. η; ℓ::Integer is not differentiable) ────────────────

function w(ℓ::Integer, η_::Dual{T}) where {T}
    ηv = value(η_)
    Dual{T}(w(ℓ, ηv), _fd_deriv(x -> w(ℓ, x), ηv) * partials(η_))
end

# ── Fresnel (analytic derivatives) ─────────────────────────────────────────────

function FresnelC(z::Dual{T}) where {T}
    zv = value(z)
    Dual{T}(FresnelC(zv), cos((π / 2) * zv^2) * partials(z))
end

function FresnelS(z::Dual{T}) where {T}
    zv = value(z)
    Dual{T}(FresnelS(zv), sin((π / 2) * zv^2) * partials(z))
end

function FresnelE(z::Dual{T}) where {T}
    zv = value(z)
    y = FresnelE(zv)
    dy = exp(im * (π / 2) * zv^2)
    complex(Dual{T}(real(y), real(dy) * partials(z)), Dual{T}(imag(y), imag(dy) * partials(z)))
end

fresnel(z::Dual{T}) where {T} = (FresnelC(z), FresnelS(z), FresnelE(z))

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
    Dual{T}(y, dθ * partials(θv))
end

# ── Fermi-Dirac (FD) ───────────────────────────────────────────────────────────

function FermiDiracIntegral(j::Dual{T}, x::Real) where {T}
    jv = value(j)
    Dual{T}(FermiDiracIntegral(jv, x), _fd_deriv(t -> FermiDiracIntegral(t, x), jv) * partials(j))
end

function FermiDiracIntegral(j::Real, x::Dual{T}) where {T}
    xv = value(x)
    Dual{T}(FermiDiracIntegral(j, xv), _fd_deriv(t -> FermiDiracIntegral(j, t), xv) * partials(x))
end

function FermiDiracIntegral(j::Dual{T}, x::Dual{T}) where {T}
    jv, xv = value(j), value(x)
    y = FermiDiracIntegral(jv, xv)
    Dual{T}(
        y,
        _fd_deriv(t -> FermiDiracIntegral(t, xv), jv) * partials(j) +
            _fd_deriv(t -> FermiDiracIntegral(jv, t), xv) * partials(x),
    )
end

function FermiDiracIntegralNorm(j::Dual{T}, x::Real) where {T}
    jv = value(j)
    Dual{T}(FermiDiracIntegralNorm(jv, x), _fd_deriv(t -> FermiDiracIntegralNorm(t, x), jv) * partials(j))
end

function FermiDiracIntegralNorm(j::Real, x::Dual{T}) where {T}
    xv = value(x)
    Dual{T}(FermiDiracIntegralNorm(j, xv), _fd_deriv(t -> FermiDiracIntegralNorm(j, t), xv) * partials(x))
end

function FermiDiracIntegralNorm(j::Dual{T}, x::Dual{T}) where {T}
    jv, xv = value(j), value(x)
    y = FermiDiracIntegralNorm(jv, xv)
    Dual{T}(
        y,
        _fd_deriv(t -> FermiDiracIntegralNorm(t, xv), jv) * partials(j) +
            _fd_deriv(t -> FermiDiracIntegralNorm(jv, t), xv) * partials(x),
    )
end

# ── MarcumQ (analytic dQ/db = dQdb; FD for M and a) ───────────────────────────

function MarcumQ(M::Dual{T}, a::Real, b::Real) where {T}
    Mv = value(M)
    Dual{T}(MarcumQ(Mv, a, b), _fd_deriv(x -> MarcumQ(x, a, b), Mv) * partials(M))
end

function MarcumQ(M::Real, a::Dual{T}, b::Real) where {T}
    av = value(a)
    Dual{T}(MarcumQ(M, av, b), _fd_deriv(x -> MarcumQ(M, x, b), av) * partials(a))
end

function MarcumQ(M::Real, a::Real, b::Dual{T}) where {T}
    bv = value(b)
    Dual{T}(MarcumQ(M, a, bv), dQdb(M, a, bv) * partials(b))
end

function MarcumQ(M::Dual{T}, a::Dual{T}, b::Real) where {T}
    Mv, av = value(M), value(a)
    y = MarcumQ(Mv, av, b)
    Dual{T}(
        y,
        _fd_deriv(x -> MarcumQ(x, av, b), Mv) * partials(M) +
            _fd_deriv(x -> MarcumQ(Mv, x, b), av) * partials(a),
    )
end

function MarcumQ(M::Dual{T}, a::Real, b::Dual{T}) where {T}
    Mv, bv = value(M), value(b)
    y = MarcumQ(Mv, a, bv)
    Dual{T}(y, _fd_deriv(x -> MarcumQ(x, a, bv), Mv) * partials(M) + dQdb(Mv, a, bv) * partials(b))
end

function MarcumQ(M::Real, a::Dual{T}, b::Dual{T}) where {T}
    av, bv = value(a), value(b)
    y = MarcumQ(M, av, bv)
    Dual{T}(y, _fd_deriv(x -> MarcumQ(M, x, bv), av) * partials(a) + dQdb(M, av, bv) * partials(b))
end

function MarcumQ(M::Dual{T}, a::Dual{T}, b::Dual{T}) where {T}
    Mv, av, bv = value(M), value(a), value(b)
    y = MarcumQ(Mv, av, bv)
    Dual{T}(
        y,
        _fd_deriv(x -> MarcumQ(x, av, bv), Mv) * partials(M) +
            _fd_deriv(x -> MarcumQ(Mv, x, bv), av) * partials(a) +
            dQdb(Mv, av, bv) * partials(b),
    )
end

# MarcumQ(a, b) — 2-arg convenience form (M=1 fixed)

function MarcumQ(a::Dual{T}, b::Real) where {T}
    av = value(a)
    Dual{T}(MarcumQ(av, b), _fd_deriv(x -> MarcumQ(x, b), av) * partials(a))
end

function MarcumQ(a::Real, b::Dual{T}) where {T}
    bv = value(b)
    Dual{T}(MarcumQ(a, bv), dQdb(a, bv) * partials(b))
end

function MarcumQ(a::Dual{T}, b::Dual{T}) where {T}
    av, bv = value(a), value(b)
    y = MarcumQ(av, bv)
    Dual{T}(y, _fd_deriv(x -> MarcumQ(x, bv), av) * partials(a) + dQdb(av, bv) * partials(b))
end

# ── dQdb (FD for all arguments) ────────────────────────────────────────────────

function dQdb(M::Dual{T}, a::Real, b::Real) where {T}
    Mv = value(M)
    Dual{T}(dQdb(Mv, a, b), _fd_deriv(x -> dQdb(x, a, b), Mv) * partials(M))
end

function dQdb(M::Real, a::Dual{T}, b::Real) where {T}
    av = value(a)
    Dual{T}(dQdb(M, av, b), _fd_deriv(x -> dQdb(M, x, b), av) * partials(a))
end

function dQdb(M::Real, a::Real, b::Dual{T}) where {T}
    bv = value(b)
    Dual{T}(dQdb(M, a, bv), _fd_deriv(x -> dQdb(M, a, x), bv) * partials(b))
end

function dQdb(M::Dual{T}, a::Dual{T}, b::Real) where {T}
    Mv, av = value(M), value(a)
    y = dQdb(Mv, av, b)
    Dual{T}(y, _fd_deriv(x -> dQdb(x, av, b), Mv) * partials(M) + _fd_deriv(x -> dQdb(Mv, x, b), av) * partials(a))
end

function dQdb(M::Dual{T}, a::Real, b::Dual{T}) where {T}
    Mv, bv = value(M), value(b)
    y = dQdb(Mv, a, bv)
    Dual{T}(y, _fd_deriv(x -> dQdb(x, a, bv), Mv) * partials(M) + _fd_deriv(x -> dQdb(Mv, a, x), bv) * partials(b))
end

function dQdb(M::Real, a::Dual{T}, b::Dual{T}) where {T}
    av, bv = value(a), value(b)
    y = dQdb(M, av, bv)
    Dual{T}(y, _fd_deriv(x -> dQdb(M, x, bv), av) * partials(a) + _fd_deriv(x -> dQdb(M, av, x), bv) * partials(b))
end

function dQdb(M::Dual{T}, a::Dual{T}, b::Dual{T}) where {T}
    Mv, av, bv = value(M), value(a), value(b)
    y = dQdb(Mv, av, bv)
    Dual{T}(
        y,
        _fd_deriv(x -> dQdb(x, av, bv), Mv) * partials(M) +
            _fd_deriv(x -> dQdb(Mv, x, bv), av) * partials(a) +
            _fd_deriv(x -> dQdb(Mv, av, x), bv) * partials(b),
    )
end

# dQdb(a, b) — 2-arg convenience form (M=1 fixed), FD

function dQdb(a::Dual{T}, b::Real) where {T}
    av = value(a)
    Dual{T}(dQdb(av, b), _fd_deriv(x -> dQdb(x, b), av) * partials(a))
end

function dQdb(a::Real, b::Dual{T}) where {T}
    bv = value(b)
    Dual{T}(dQdb(a, bv), _fd_deriv(x -> dQdb(a, x), bv) * partials(b))
end

function dQdb(a::Dual{T}, b::Dual{T}) where {T}
    av, bv = value(a), value(b)
    y = dQdb(av, bv)
    Dual{T}(y, _fd_deriv(x -> dQdb(x, bv), av) * partials(a) + _fd_deriv(x -> dQdb(av, x), bv) * partials(b))
end

# ── Debye function (FD) ────────────────────────────────────────────────────────

# 2-arg form: debye_function(β, x)

function debye_function(β::Dual{T}, x::Real; tol = 1.0e-35, max_terms = 2000) where {T}
    βv = value(β)
    y = debye_function(βv, x; tol = tol, max_terms = max_terms)
    dβ = _fd_deriv(t -> debye_function(t, x; tol = tol, max_terms = max_terms), βv)
    Dual{T}(y, dβ * partials(β))
end

function debye_function(β::Real, x::Dual{T}; tol = 1.0e-35, max_terms = 2000) where {T}
    xv = value(x)
    y = debye_function(β, xv; tol = tol, max_terms = max_terms)
    dx = _fd_deriv(t -> debye_function(β, t; tol = tol, max_terms = max_terms), xv)
    Dual{T}(y, dx * partials(x))
end

function debye_function(β::Dual{T}, x::Dual{T}; tol = 1.0e-35, max_terms = 2000) where {T}
    βv, xv = value(β), value(x)
    y = debye_function(βv, xv; tol = tol, max_terms = max_terms)
    dβ = _fd_deriv(t -> debye_function(t, xv; tol = tol, max_terms = max_terms), βv)
    dx = _fd_deriv(t -> debye_function(βv, t; tol = tol, max_terms = max_terms), xv)
    Dual{T}(y, dβ * partials(β) + dx * partials(x))
end

# 3-arg form: debye_function(n, β, x)

function debye_function(n::Dual{T}, β::Real, x::Real; tol = 1.0e-35, max_terms = 2000) where {T}
    nv = value(n)
    y = debye_function(nv, β, x; tol = tol, max_terms = max_terms)
    dn = _fd_deriv(t -> debye_function(t, β, x; tol = tol, max_terms = max_terms), nv)
    Dual{T}(y, dn * partials(n))
end

function debye_function(n::Real, β::Dual{T}, x::Real; tol = 1.0e-35, max_terms = 2000) where {T}
    βv = value(β)
    y = debye_function(n, βv, x; tol = tol, max_terms = max_terms)
    dβ = _fd_deriv(t -> debye_function(n, t, x; tol = tol, max_terms = max_terms), βv)
    Dual{T}(y, dβ * partials(β))
end

function debye_function(n::Real, β::Real, x::Dual{T}; tol = 1.0e-35, max_terms = 2000) where {T}
    xv = value(x)
    y = debye_function(n, β, xv; tol = tol, max_terms = max_terms)
    dx = _fd_deriv(t -> debye_function(n, β, t; tol = tol, max_terms = max_terms), xv)
    Dual{T}(y, dx * partials(x))
end

function debye_function(n::Dual{T}, β::Dual{T}, x::Real; tol = 1.0e-35, max_terms = 2000) where {T}
    nv, βv = value(n), value(β)
    y = debye_function(nv, βv, x; tol = tol, max_terms = max_terms)
    dn = _fd_deriv(t -> debye_function(t, βv, x; tol = tol, max_terms = max_terms), nv)
    dβ = _fd_deriv(t -> debye_function(nv, t, x; tol = tol, max_terms = max_terms), βv)
    Dual{T}(y, dn * partials(n) + dβ * partials(β))
end

function debye_function(n::Dual{T}, β::Real, x::Dual{T}; tol = 1.0e-35, max_terms = 2000) where {T}
    nv, xv = value(n), value(x)
    y = debye_function(nv, β, xv; tol = tol, max_terms = max_terms)
    dn = _fd_deriv(t -> debye_function(t, β, xv; tol = tol, max_terms = max_terms), nv)
    dx = _fd_deriv(t -> debye_function(nv, β, t; tol = tol, max_terms = max_terms), xv)
    Dual{T}(y, dn * partials(n) + dx * partials(x))
end

function debye_function(n::Real, β::Dual{T}, x::Dual{T}; tol = 1.0e-35, max_terms = 2000) where {T}
    βv, xv = value(β), value(x)
    y = debye_function(n, βv, xv; tol = tol, max_terms = max_terms)
    dβ = _fd_deriv(t -> debye_function(n, t, xv; tol = tol, max_terms = max_terms), βv)
    dx = _fd_deriv(t -> debye_function(n, βv, t; tol = tol, max_terms = max_terms), xv)
    Dual{T}(y, dβ * partials(β) + dx * partials(x))
end

function debye_function(n::Dual{T}, β::Dual{T}, x::Dual{T}; tol = 1.0e-35, max_terms = 2000) where {T}
    nv, βv, xv = value(n), value(β), value(x)
    y = debye_function(nv, βv, xv; tol = tol, max_terms = max_terms)
    dn = _fd_deriv(t -> debye_function(t, βv, xv; tol = tol, max_terms = max_terms), nv)
    dβ = _fd_deriv(t -> debye_function(nv, t, xv; tol = tol, max_terms = max_terms), βv)
    dx = _fd_deriv(t -> debye_function(nv, βv, t; tol = tol, max_terms = max_terms), xv)
    Dual{T}(y, dn * partials(n) + dβ * partials(β) + dx * partials(x))
end

# ── Parabolic cylinder: U, V, W ────────────────────────────────────────────────
# ∂U/∂x = dU(a, x),  ∂V/∂x = dV(a, x),  ∂W/∂x = dW(a, x)

function U(a::Dual{T}, x::Real) where {T}
    av = value(a)
    Dual{T}(U(av, x), _fd_deriv(t -> U(t, x), av) * partials(a))
end

function U(a::Real, x::Dual{T}) where {T}
    xv = value(x)
    Dual{T}(U(a, xv), dU(a, xv) * partials(x))
end

function U(a::Dual{T}, x::Dual{T}) where {T}
    av, xv = value(a), value(x)
    Dual{T}(U(av, xv), _fd_deriv(t -> U(t, xv), av) * partials(a) + dU(av, xv) * partials(x))
end

function V(a::Dual{T}, x::Real) where {T}
    av = value(a)
    Dual{T}(V(av, x), _fd_deriv(t -> V(t, x), av) * partials(a))
end

function V(a::Real, x::Dual{T}) where {T}
    xv = value(x)
    Dual{T}(V(a, xv), dV(a, xv) * partials(x))
end

function V(a::Dual{T}, x::Dual{T}) where {T}
    av, xv = value(a), value(x)
    Dual{T}(V(av, xv), _fd_deriv(t -> V(t, xv), av) * partials(a) + dV(av, xv) * partials(x))
end

function W(a::Dual{T}, x::Real) where {T}
    av = value(a)
    Dual{T}(W(av, x), _fd_deriv(t -> W(t, x), av) * partials(a))
end

function W(a::Real, x::Dual{T}) where {T}
    xv = value(x)
    Dual{T}(W(a, xv), dW(a, xv) * partials(x))
end

function W(a::Dual{T}, x::Dual{T}) where {T}
    av, xv = value(a), value(x)
    Dual{T}(W(av, xv), _fd_deriv(t -> W(t, xv), av) * partials(a) + dW(av, xv) * partials(x))
end

# ── Parabolic cylinder: dU, dV ─────────────────────────────────────────────────
# d(dU)/dx = (x²/4 + a) * U(a, x)  from the parabolic cylinder ODE U'' = (x²/4 + a)U

function dU(a::Dual{T}, x::Real) where {T}
    av = value(a)
    Dual{T}(dU(av, x), _fd_deriv(t -> dU(t, x), av) * partials(a))
end

function dU(a::Real, x::Dual{T}) where {T}
    xv = value(x)
    Dual{T}(dU(a, xv), (xv^2 / 4 + a) * U(a, xv) * partials(x))
end

function dU(a::Dual{T}, x::Dual{T}) where {T}
    av, xv = value(a), value(x)
    Dual{T}(dU(av, xv), _fd_deriv(t -> dU(t, xv), av) * partials(a) + (xv^2 / 4 + av) * U(av, xv) * partials(x))
end

function dV(a::Dual{T}, x::Real) where {T}
    av = value(a)
    Dual{T}(dV(av, x), _fd_deriv(t -> dV(t, x), av) * partials(a))
end

function dV(a::Real, x::Dual{T}) where {T}
    xv = value(x)
    Dual{T}(dV(a, xv), (xv^2 / 4 + a) * V(a, xv) * partials(x))
end

function dV(a::Dual{T}, x::Dual{T}) where {T}
    av, xv = value(a), value(x)
    Dual{T}(dV(av, xv), _fd_deriv(t -> dV(t, xv), av) * partials(a) + (xv^2 / 4 + av) * V(av, xv) * partials(x))
end

# ── Parabolic cylinder: dW (FD for both arguments) ────────────────────────────

function dW(a::Dual{T}, x::Real) where {T}
    av = value(a)
    Dual{T}(dW(av, x), _fd_deriv(t -> dW(t, x), av) * partials(a))
end

function dW(a::Real, x::Dual{T}) where {T}
    xv = value(x)
    Dual{T}(dW(a, xv), _fd_deriv(t -> dW(a, t), xv) * partials(x))
end

function dW(a::Dual{T}, x::Dual{T}) where {T}
    av, xv = value(a), value(x)
    Dual{T}(dW(av, xv), _fd_deriv(t -> dW(t, xv), av) * partials(a) + _fd_deriv(t -> dW(av, t), xv) * partials(x))
end

end
