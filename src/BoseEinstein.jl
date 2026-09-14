using SpecialFunctions: gamma, loggamma, zeta

@doc raw"""
    BoseEinsteinIntegralNorm(k, η)

Compute the normalized Bose–Einstein integral

```math
B_k(\eta) = \frac{1}{\Gamma(k+1)}\int_0^\infty
\frac{t^k}{\exp(t-\eta)-1}\,dt = \operatorname{Li}_{k+1}(e^\eta).
```

Accepts integer and half-integer orders `k ≥ -9/2` and real `η ≤ 0`.
For `k ≤ -1`, the polylogarithm defines the continuation by differentiation;
the displayed integral itself requires `k > -1`. At `η = 0`, returns
`zeta(k+1)` for `k > 0` and `Inf` otherwise. Returns zero at `η = -Inf`.
Invalid orders, positive `η`, and NaN inputs throw `DomainError`.

Uses Fukushima's minimax rational approximations for half-integer orders
`-9/2:1:39/2` and integer orders `1:19`, elementary formulas for nonpositive
integers, and a convergent series for higher orders. `Float32` and `Float64`
use the double precision coefficients; `BigFloat` uses convergent series at
the working precision. The ForwardDiff
extension supports differentiation with respect to `η`, but not the discrete
order `k`.

# Examples
```jldoctest
julia> round(BoseEinsteinIntegralNorm(0.5, -1.0); digits=12)
0.4284407346
```

# Reference
T. Fukushima (2020 preprint), *Analytical computation of Bose–Einstein integral of half
integer orders, −9/2, −7/2, ⋯, and 39/2, and integer orders, 1, 2, ⋯, and 19,
by minimax rational function approximations*, Eqs. (11)–(24), Tables A.3–A.48,
[doi:10.13140/RG.2.2.21720.65283](https://doi.org/10.13140/RG.2.2.21720.65283).

See also [`BoseEinsteinIntegral`](@ref), [`FermiDiracIntegralNorm`](@ref).
"""
function BoseEinsteinIntegralNorm(k::Real, η::Real)
    T = float(promote_type(typeof(k), typeof(η)))
    kT, ηT = T(k), T(η)
    _bose_check_domain(kT, ηT)
    if T === Float16 || T === Float32
        return T(_bose_normalized(Float64(kT), Float64(ηT)))
    end
    return _bose_normalized(kT, ηT)
end

@doc raw"""
    BoseEinsteinIntegral(k, η)

Compute the unnormalized Bose–Einstein integral

```math
\int_0^\infty \frac{t^k}{\exp(t-\eta)-1}\,dt
= \Gamma(k+1) B_k(\eta).
```

Accepts integer and half-integer orders `k > -1` and real `η ≤ 0`.
At `η = 0`, returns `gamma(k+1)*zeta(k+1)` for `k > 0`, and `Inf`
for `k = -1/2` or `k = 0`. Returns zero at `η = -Inf`. Invalid inputs
throw `DomainError`. Supports `Float32`, `Float64`, `BigFloat`, dot
broadcasting, and ForwardDiff differentiation with respect to `η`.

The naming follows [`FermiDiracIntegral`](@ref). For normalized values and
continuation to negative orders, use [`BoseEinsteinIntegralNorm`](@ref).
"""
function BoseEinsteinIntegral(k::Real, η::Real)
    T = float(promote_type(typeof(k), typeof(η)))
    kT, ηT = T(k), T(η)
    _bose_check_domain(kT, ηT)
    kT > -1 || throw(DomainError(k, "the unnormalized integral requires k > -1"))
    if T === Float16 || T === Float32
        return T(_bose_unnormalized(Float64(kT), Float64(ηT)))
    end
    return _bose_unnormalized(kT, ηT)
end

function _bose_check_domain(k, η)
    isfinite(k) && k >= -4.5 && (isinteger(k) || isinteger(k - 1 / 2)) ||
        throw(DomainError(k, "k must be an integer or half-integer at least -9/2"))
    η <= 0 || throw(DomainError(η, "η must be nonpositive"))
    return nothing
end

function _bose_normalized(k::T, η::T) where {T <: AbstractFloat}
    η == -Inf && return zero(T)
    η == 0 && return k > 0 ? T(zeta(k + one(T))) : T(Inf)
    if k == 0
        return η < -log(T(2)) ? -log1p(-exp(η)) : -log(-expm1(η))
    elseif isinteger(k) && -4 <= k < 0
        z, q = exp(η), -expm1(η)
        k == -1 && return z / q
        k == -2 && return z / q^2
        k == -3 && return z * (one(T) + z) / q^3
        return z * (one(T) + 4z + z^2) / q^4
    elseif T === Float64 && -4.5 <= k < 20
        return _bose_minimax(k, η)
    end
    return _bose_series(k, η)
end

function _bose_unnormalized(k::T, η::T, order::T = k) where {T <: AbstractFloat}
    η == -Inf && return zero(T)
    b = _bose_normalized(order, η)
    # Avoid an overflowing gamma or an underflowing exp(η) before their product.
    if (T === BigFloat || k <= 170.5) && b >= floatmin(T)
        scale = gamma(k + one(T))
        isfinite(scale) && return scale * b
    end
    logb = η < -1 ? η + log(_bose_fugacity_scaled(order, η)) : log(b)
    return exp(loggamma(k + one(T)) + logb)
end

# Eq. (11), divided by exp(η) so the first term is exactly one. For η near
# zero and k > 0, the integral bound on the p-series controls the positive tail.
function _bose_fugacity_scaled(k::T, η::T) where {T <: AbstractFloat}
    z = exp(η)
    q = -expm1(η)
    power, total, correction = one(T), one(T), zero(T)
    for n in 2:100000
        power *= z
        term = power * T(n)^(-k - one(T))
        adjusted = term - correction
        updated = total + adjusted
        correction = (updated - total) - adjusted
        total = updated
        tail = if k >= -1
            geometric = term * z / q
            k > 0 ? min(geometric, term * T(n) / k) : geometric
        else
            ratio = z * (one(T) + inv(T(n)))^(-k - one(T))
            ratio < 1 ? term * ratio / (one(T) - ratio) : T(Inf)
        end
        tail <= eps(T) * total / 4 && return total
    end
    error("Bose–Einstein fugacity series did not converge")
end

function _bose_series(k::T, η::T) where {T <: AbstractFloat}
    η == -Inf && return zero(T)
    η == 0 && return k > 0 ? T(zeta(k + one(T))) : T(Inf)
    # High orders converge rapidly even at exp(η) ≈ 1.
    if η <= -1 || k >= precision(T)
        return exp(η) * _bose_fugacity_scaled(k, η)
    end
    # Float64 orders ≥ 20 need only a few terms of Eq. (11).
    T === Float64 && k >= 20 && return exp(η) * _bose_fugacity_scaled(k, η)

    integer_order = isinteger(k) && k >= 0
    total = if integer_order
        n = Int(k)
        harmonic = sum(inv(T(j)) for j in 1:n; init = zero(T))
        # The sign follows η^k, as in Eq. (23); Eq. (14) prints (-η)^k.
        η^n / gamma(k + one(T)) * (harmonic - log(-η))
    else
        gamma(-k) * (-η)^k
    end
    isinf(total) && return total
    coefficient, correction = one(T), zero(T)
    small_terms = 0
    for m in 0:100000
        if !integer_order || m != k
            term = zeta(k + one(T) - m) * coefficient
            adjusted = term - correction
            updated = total + adjusted
            correction = (updated - total) - adjusted
            total = updated
            small_terms = abs(term) <= eps(T) * abs(total) / 4 ? small_terms + 1 : 0
            # ζ vanishes at negative even integers; never stop on one zero term.
            m > k + 1 && small_terms >= 3 && return total
        end
        coefficient *= η / (m + 1)
    end
    error("Bose–Einstein expansion at η = 0 did not converge")
end

include("BoseEinsteinCoefficients.jl")
