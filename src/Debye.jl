using QuadGK: quadgk

export debye_function

"""
    debye_function(n::T, β::T, x::T; tol=1e-35, max_terms=2000) where {T <: AbstractFloat}

Compute `n/x^n * ∫₀ˣ t^n/(exp(t)-1)^β dt` for finite `n > 0`, `0 < β < n+1`,
and `x ≥ 0`. At `x = 0`, returns the limiting value: `0` for `β < 1`, `1` for
`β = 1`, and `Inf` for `β > 1`. Returns `0` at `x = Inf`. Supports any `AbstractFloat` type
(e.g., `Float32`, `Float64`, `BigFloat`).
Array broadcasting is supported: any one argument may be an `AbstractArray`.

Throws `ArgumentError` for invalid parameters, tolerances, or work limits.
Throws `ErrorException` if quadrature cannot meet the requested tolerance within the work limit.

- `tol`: positive finite relative quadrature tolerance (default `1e-35`), floored at `8eps(T)`
- `max_terms`: positive integer limit on quadrature subintervals (default `2000`)

References:
- [Debye function](https://en.wikipedia.org/wiki/Debye_function)
- [Paper](https://doi.org/10.1007/s10765-007-0256-1)
"""
function debye_function(n::T, β::T, x::T; tol = 1.0e-35, max_terms = 2000) where {T <: AbstractFloat}
    isfinite(n) && n > 0 || throw(ArgumentError("n must be positive and finite, got n = $n"))
    q = β <= one(T) ? n + (one(T) - β) : (n - β) + one(T)
    isfinite(β) && β > 0 && q > 0 || throw(ArgumentError("β must satisfy 0 < β < n+1, got β = $β"))
    x >= 0 || throw(ArgumentError("x must be non-negative, got x = $x"))
    isfinite(tol) && tol > 0 || throw(ArgumentError("tol must be positive and finite, got tol = $tol"))
    max_terms isa Integer && max_terms > 0 || throw(ArgumentError("max_terms must be a positive integer"))
    x == 0 && return β < one(T) ? zero(T) : β == one(T) ? one(T) : T(Inf)
    isinf(x) && return zero(T)

    m = min(one(T), q)
    scale = min(x, max(one(T), n / β))
    # t = scale*s/(1-s) keeps even a very large upper limit on a bounded interval.
    # v = s^m removes the integrable singularity at the origin when q < 1.
    integrand = function (v)
        s = v^(one(T) / m)
        ratio = s / (one(T) - s)
        t = scale * ratio
        isinf(t) && return zero(T)
        log_bose = if t == 0
            zero(T)
        elseif t <= one(T)
            log(t / expm1(t))
        else
            log(t) - t - log(-expm1(-t))
        end
        log_power = q == m ? zero(T) : (q - m) * log(ratio)
        return exp(log_power - (m + one(T)) * log1p(-s) + β * log_bose)
    end
    rtol = max(T(tol), 8eps(T))
    integral, error = quadgk(
        integrand, zero(T), (one(T) / (one(T) + scale / x))^m;
        rtol, maxevals = 15 * (2big(max_terms) - 1)
    )
    isfinite(integral) && integral > 0 && error <= rtol * integral ||
        throw(ErrorException("Debye quadrature did not converge within max_terms = $max_terms"))
    result = (n / m) * integral * scale^(one(T) - β) * (scale / x)^n
    isfinite(result) && result > 0 && return result
    # Recover representable results when separate scaling factors overflow or underflow.
    return exp(log(n) - log(m) + (one(T) - β) * log(scale) + n * (log(scale) - log(x)) + log(integral))
end

function debye_function(n::Real, β::Real, x::Real; tol = 1.0e-35, max_terms = 2000)
    T = float(promote_type(typeof(n), typeof(β), typeof(x)))
    return debye_function(T(n), T(β), T(x); tol = tol, max_terms = max_terms)
end

function debye_function(n::Real, β::Real, x::AbstractArray{<:Real}; tol = 1.0e-35, max_terms = 2000)
    return [debye_function(n, β, xᵢ; tol = tol, max_terms = max_terms) for xᵢ in x]
end

function debye_function(n::Real, β::AbstractArray{<:Real}, x::Real; tol = 1.0e-35, max_terms = 2000)
    return [debye_function(n, βᵢ, x; tol = tol, max_terms = max_terms) for βᵢ in β]
end

function debye_function(n::AbstractArray{<:Real}, β::Real, x::Real; tol = 1.0e-35, max_terms = 2000)
    return [debye_function(nᵢ, β, x; tol = tol, max_terms = max_terms) for nᵢ in n]
end

function debye_function(β::Real, x::Real; tol = 1.0e-35, max_terms = 2000)
    return debye_function(one(float(promote_type(typeof(β), typeof(x)))), β, x; tol = tol, max_terms = max_terms)
end
