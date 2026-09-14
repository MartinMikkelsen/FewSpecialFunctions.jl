using SpecialFunctions: gamma, loggamma

export U, V, W, dU, dV, dW

# Older supported Julia versions have process-global MPFR precision.
const _cylinder_precision_lock = ReentrantLock()

# Values at zero, including reciprocal-gamma zeros (DLMF 12.2 and 12.14).
_cylinder_rgamma(x) = x <= 0 && isinteger(x) ? zero(x) : inv(gamma(x))

function _cylinder_initial(a::T, kind::Symbol) where {T <: AbstractFloat}
    if kind === :U
        c = sqrt(T(π)) * T(2)^(-a / 2 - T(0.25))
        return c * _cylinder_rgamma(a / 2 + T(0.75)),
            -sqrt(T(2)) * c * _cylinder_rgamma(a / 2 + T(0.25))
    elseif kind === :V
        return T(2)^(a / 2 + T(0.25)) * sinpi(a / 2 + T(0.25)) * _cylinder_rgamma(T(0.75) - a / 2),
            T(2)^(a / 2 + T(0.75)) * sinpi(a / 2 + T(0.75)) * _cylinder_rgamma(T(0.25) - a / 2)
    else
        l = real(loggamma(Complex{T}(T(0.25), a / 2)) - loggamma(Complex{T}(T(0.75), a / 2))) / 2
        return T(2)^(-T(0.75)) * exp(l), -T(2)^(-T(0.25)) * exp(-l)
    end
end

# Sum actual Taylor terms of y'' = (a ± x²/4)y, together with their
# derivatives. Keeping the factorial in the recurrence avoids overflowing
# coefficient arrays, even when the final terms are small.
function _cylinder_series(a::T, x::T, kind::Symbol) where {T <: AbstractFloat}
    y0, y1 = _cylinder_initial(a, kind)
    iszero(x) && return y0, y1
    ep, e, op, o = zero(T), y0, zero(T), x * y1
    y, dy = e + o, y1
    A, B = a * x^2, (kind === :W ? -one(T) : one(T)) * x^4 / 4
    small = 0
    for k in 1:10000
        ep, e = e, (A * e + B * ep) / (T(2k) * T(2k - 1))
        op, o = o, (A * o + B * op) / (T(2k + 1) * T(2k))
        dyterm = (T(2k) * e + T(2k + 1) * o) / x
        y += e + o
        dy += dyterm
        value_small = abs(e) + abs(o) <= eps(T) * abs(y)
        derivative_small = (T(2k) * abs(e) + T(2k + 1) * abs(o)) / abs(x) <= eps(T) * abs(dy)
        small = value_small && derivative_small ? small + 1 : 0
        small >= 3 && return y, dy
    end
    throw(ErrorException("Parabolic cylinder series did not converge"))
end

# Growing and decaying solutions can cancel in a Taylor sum. Reserve bits
# for both exponential scales, exp(±(x²/4 + sqrt(abs(a))*abs(x))), then round
# back to the caller's precision. This fallback trades speed for accuracy.
function _cylinder_series_guarded(a::T, x::T, kind::Symbol, scaled::Bool = false) where {T <: AbstractFloat}
    isfinite(a) && isfinite(x) || throw(DomainError((a, x), "cylinder parameters must be finite"))
    iszero(x) && !scaled && return _cylinder_initial(a, kind)
    p = precision(x)
    extra = ceil(Int, (abs(x)^2 / 2 + 2sqrt(abs(a)) * abs(x)) / log(T(2))) + 32
    values = lock(_cylinder_precision_lock) do
        setprecision(BigFloat, p + extra) do
            ab, xb = BigFloat(a), BigFloat(x)
            result = _cylinder_series(ab, xb, kind)
            if scaled
                logscale = _cylinder_logscale(ab, xb) * (kind === :U ? 1 : -1)
                return map(y -> _cylinder_scaled(logscale, y), result)
            end
            return result
        end
    end
    return map(v -> T === BigFloat ? BigFloat(v; precision = p) : T(v), values)
end

# DLMF 12.9.1. Use an asymptotic expansion only if its decreasing terms
# reach the target precision; otherwise return nothing and use the series.
function _cylinder_scaled(logscale::T, y::T) where {T <: AbstractFloat}
    scale = exp(logscale)
    if scale < floatmin(T) || !isfinite(scale)
        return iszero(y) ? y : copysign(exp(logscale + log(abs(y))), y)
    end
    return scale * y
end

function _cylinder_u_asymptotic(a::T, x::T, scaled::Bool = false) where {T <: AbstractFloat}
    r, s, ds = one(T), one(T), zero(T)
    for k in 1:1000
        next = -r * (a + T(2k) - T(1.5)) * (a + T(2k) - T(0.5)) / (T(2k) * x^2)
        abs(next) > abs(r) && return nothing
        s += next
        ds -= T(2k) * next / x
        r = next
        if abs(r) <= eps(T) * abs(s) / 16
            logscale = scaled ? _cylinder_asymptotic_scale(a, x, :U) : -x^2 / 4 - (a + T(0.5)) * log(x)
            return _cylinder_scaled(logscale, s), _cylinder_scaled(logscale, ds - (x / 2 + (a + T(0.5)) / x) * s)
        end
    end
    return nothing
end

# Gil, Segura & Temme (2006), equations (7), (11)–(13). The absolute
# value of the saddle-point factor is constant in the oscillatory region.
function _cylinder_logscale(a::T, x::T) where {T <: AbstractFloat}
    iszero(a) && return x^2 / 4
    q = x^2 / 4 + a
    q <= 0 && return a * (log(abs(a)) - 1) / 2
    d = sqrt(q)
    return a * log(x / 2 + d) + x * d / 2 - a / 2
end

# Combine the saddle-point and Poincaré factors without subtracting x²/4.
function _cylinder_asymptotic_scale(a::T, x::T, kind::Symbol) where {T <: AbstractFloat}
    q = a / (sqrt(x^2 / 4 + a) + x / 2)
    correction = a * log1p(q / x) - q^2 / 2
    return (kind === :U ? correction : -correction) - log(x) / 2
end

function _cylinder_v_asymptotic_scaled(a::T, x::T) where {T <: AbstractFloat}
    r, s, ds = one(T), one(T), zero(T)
    for k in 1:1000
        next = r * (T(2k) - T(1.5) - a) * (T(2k) - T(0.5) - a) / (T(2k) * x^2)
        abs(next) > abs(r) && return nothing
        s += next
        ds -= T(2k) * next / x
        r = next
        if abs(r) <= eps(T) * abs(s) / 16
            logscale = _cylinder_asymptotic_scale(a, x, :V) + log(T(2) / T(π)) / 2
            return _cylinder_scaled(logscale, s), _cylinder_scaled(logscale, ds + (x / 2 + (a - T(0.5)) / x) * s)
        end
    end
    return nothing
end

function _cylinder_scaled_pair(a::T, x::T, kind::Symbol) where {T <: AbstractFloat}
    isfinite(a) && isfinite(x) && x >= 0 || throw(DomainError((a, x), "scaled cylinder functions require finite a and x ≥ 0"))
    if x > 1 && x^2 / 4 + a > 0
        result = kind === :U ? _cylinder_u_asymptotic(a, x, true) : _cylinder_v_asymptotic_scaled(a, x)
        result === nothing || return result
    end
    return _cylinder_series_guarded(a, x, kind, true)
end

"""
    U_scaled(a::Real, x::Real)

Compute `F(a,x) * U(a,x)` for finite real `a` and `x ≥ 0`, evaluating the
scaling before rounding to the output type. Here `F = exp(L)`, with
`L = a*log(x/2 + sqrt(x²/4+a)) + x*sqrt(x²/4+a)/2 - a/2` in the
nonoscillatory region, `L = a*(log(abs(a))-1)/2` in the oscillatory region,
and `L = x²/4` when `a = 0`.

Supports Float32, Float64 and BigFloat. Invalid arguments raise `DomainError`.
The series fallback uses extra precision and can be expensive for large orders.
This scaling is from Gil, Segura & Temme (2006), equations (7), (11)–(13):
[paper](https://ir.cwi.nl/pub/14654/14654D.pdf). It removes growth/decay in both
order and argument; it is not simply multiplication by `exp(x²/4)`.
"""
U_scaled(a::T, x::T) where {T <: AbstractFloat} = first(_cylinder_scaled_pair(a, x, :U))

"""
    V_scaled(a::Real, x::Real)

Compute `V(a,x) / F(a,x)` for finite real `a` and `x ≥ 0`, with the scaling
factor defined in [`U_scaled`](@ref). Scaling is applied before rounding,
so values can remain finite when `V` overflows. Supports Float32, Float64
and BigFloat; invalid arguments raise `DomainError`.
"""
V_scaled(a::T, x::T) where {T <: AbstractFloat} = first(_cylinder_scaled_pair(a, x, :V))

"""
    ParabolicCylinderD(ν::Real, x::Real)

Compute `Dν(x) = U(-ν-1/2, x)` for finite real order and argument.
Supports Float32, Float64 and BigFloat, with mixed inputs promoted.
See [DLMF 12.2.5](https://dlmf.nist.gov/12.2.E5).
"""
function ParabolicCylinderD(ν::Real, x::Real)
    isfinite(ν) && isfinite(x) || throw(DomainError((ν, x), "cylinder parameters must be finite"))
    νf, xf = promote(float(ν), float(x))
    return U(-νf - one(νf) / 2, xf)
end

"""
    dParabolicCylinderD(ν::Real, x::Real)

Compute the argument derivative `dDν(x)/dx = dU(-ν-1/2, x)`.
The domain and type behavior match [`ParabolicCylinderD`](@ref).
"""
function dParabolicCylinderD(ν::Real, x::Real)
    isfinite(ν) && isfinite(x) || throw(DomainError((ν, x), "cylinder parameters must be finite"))
    νf, xf = promote(float(ν), float(x))
    return dU(-νf - one(νf) / 2, xf)
end

"""
    ParabolicCylinderD_scaled(ν::Real, x::Real)

Compute `U_scaled(-ν-1/2, x)` for finite real `ν` and `x ≥ 0`.
Uses the order-dependent scaling of [`U_scaled`](@ref), and supports
Float32, Float64 and BigFloat. Invalid arguments raise `DomainError`.
"""
function ParabolicCylinderD_scaled(ν::Real, x::Real)
    νf, xf = promote(float(ν), float(x))
    return U_scaled(-νf - one(νf) / 2, xf)
end

# DLMF 12.14.17–22. The coefficient of x^(-2k) is
# (-i)^k (1/2+ia)_(2k)/(k! 2^k). Differentiating the same sum keeps
# W and dW consistent; log(k) avoids subtracting nearly equal amplitudes.
function _cylinder_w_asymptotic(a::T, x::T) where {T <: AbstractFloat}
    t = abs(x)
    r = s = one(Complex{T})
    ds = zero(Complex{T})
    for k in 1:1000
        next = -im * r * Complex{T}(T(2k) - T(1.5), a) *
            Complex{T}(T(2k) - T(0.5), a) / (T(2k) * t^2)
        abs(next) > abs(r) && return nothing
        s += next
        ds -= T(2k) * next / t
        r = next
        if abs(r) <= eps(T) * abs(s) / 16
            log_inv_k = a >= 0 ? T(π) * a + log1p(sqrt(one(T) + exp(-2T(π) * a))) : asinh(exp(T(π) * a))
            phase = t^2 / 4 - a * log(t) + T(π) / 4 + imag(loggamma(Complex{T}(T(0.5), a))) / 2
            logscale = (log(T(2)) - log(t) + (x > 0 ? -log_inv_k : log_inv_k)) / 2
            oscillation = cis(phase)
            derivative = ds + Complex{T}(-one(T) / (2t), t / 2 - a / t) * s
            component = x > 0 ? real : imag
            return _cylinder_scaled(logscale, component(oscillation * s)), sign(x) * _cylinder_scaled(logscale, component(oscillation * derivative))
        end
    end
    return nothing
end

function _cylinder_u(a::T, x::T) where {T <: AbstractFloat}
    # Exact Hermite parity prevents subtracting growing solutions at negative
    # half-integer orders. Test exact equality, not an approximate order.
    if x < 0 && a <= -T(0.5) && isinteger(a + T(0.5))
        y, dy = _cylinder_u(a, -x)
        parity = -sinpi(a)
        return parity * y, -parity * dy
    end
    if x > 1
        result = _cylinder_u_asymptotic(a, x)
        result === nothing || return result
    end
    return _cylinder_series_guarded(a, x, :U)
end

function _cylinder_w(a::T, x::T) where {T <: AbstractFloat}
    if abs(x) > 1
        result = _cylinder_w_asymptotic(a, x)
        result === nothing || return result
    end
    return _cylinder_series_guarded(a, x, :W)
end

"""
    U(a::T, x::T) where {T <: AbstractFloat}

Compute the real parabolic cylinder function U(a,x). Supports Float32,
Float64, and BigFloat. Uses a precision-guarded convergent series when the
order and argument do not permit an accurate asymptotic expansion.

Reference: [DLMF, Chapter 12](https://dlmf.nist.gov/12).
"""
U(a::T, x::T) where {T <: AbstractFloat} = first(_cylinder_u(a, x))

"""
    V(a::T, x::T) where {T <: AbstractFloat}

Compute the real parabolic cylinder function V(a,x) using a convergent
series with extra working precision. Supports Float32, Float64, and BigFloat.
"""
V(a::T, x::T) where {T <: AbstractFloat} = first(_cylinder_series_guarded(a, x, :V))

"""
    W(a::T, x::T) where {T <: AbstractFloat}

Compute the real parabolic cylinder function W(a,x). Supports Float32,
Float64, and BigFloat. Uses a precision-guarded convergent series or an
asymptotic expansion whose decreasing terms reach the requested precision.

Reference: [DLMF 12.14](https://dlmf.nist.gov/12.14).
"""
W(a::T, x::T) where {T <: AbstractFloat} = first(_cylinder_w(a, x))

"""
    dU(a::T, x::T) where {T <: AbstractFloat}

Compute ∂U(a,x)/∂x from the same expansion as U. Supports Float32, Float64,
and BigFloat.
"""
dU(a::T, x::T) where {T <: AbstractFloat} = last(_cylinder_u(a, x))

"""
    dV(a::T, x::T) where {T <: AbstractFloat}

Compute ∂V(a,x)/∂x from the same convergent series as V. Supports Float32,
Float64, and BigFloat.
"""
dV(a::T, x::T) where {T <: AbstractFloat} = last(_cylinder_series_guarded(a, x, :V))

"""
    dW(a::T, x::T) where {T <: AbstractFloat}

Compute ∂W(a,x)/∂x from the same expansion as W. Supports Float32, Float64,
and BigFloat.
"""
dW(a::T, x::T) where {T <: AbstractFloat} = last(_cylinder_w(a, x))

# Extra phase bits matter for oscillatory W, particularly near its zeros.
for func in (:U, :V, :W, :dU, :dV, :dW, :U_scaled, :V_scaled)
    @eval $func(a::T, x::T) where {T <: Union{Float16, Float32}} = T($func(Float64(a), Float64(x)))
end

# Protect BigFloat asymptotics and origin values from concurrent fallbacks too.
for (func, evaluator, component) in (
        (:U, :_cylinder_u, :first), (:dU, :_cylinder_u, :last),
        (:W, :_cylinder_w, :first), (:dW, :_cylinder_w, :last),
    )
    @eval $func(a::BigFloat, x::BigFloat) = lock(() -> $component($evaluator(a, x)), _cylinder_precision_lock)
end
V(a::BigFloat, x::BigFloat) = lock(() -> first(_cylinder_series_guarded(a, x, :V)), _cylinder_precision_lock)
dV(a::BigFloat, x::BigFloat) = lock(() -> last(_cylinder_series_guarded(a, x, :V)), _cylinder_precision_lock)
U_scaled(a::BigFloat, x::BigFloat) = lock(() -> first(_cylinder_scaled_pair(a, x, :U)), _cylinder_precision_lock)
V_scaled(a::BigFloat, x::BigFloat) = lock(() -> first(_cylinder_scaled_pair(a, x, :V)), _cylinder_precision_lock)

# Promotion methods for Real inputs
for func in (:U, :V, :W, :dU, :dV, :dW, :U_scaled, :V_scaled)
    @eval function $func(a::Real, x::Real)
        T = float(promote_type(typeof(a), typeof(x)))
        return $func(T(a), T(x))
    end
end

U(a::Real, x::AbstractArray{<:Real}) = [U(a, xi) for xi in x]
U(a::AbstractArray{<:Real}, x::Real) = [U(ai, x) for ai in a]

dU(a::Real, x::AbstractArray{<:Real}) = [dU(a, xi) for xi in x]
dU(a::AbstractArray{<:Real}, x::Real) = [dU(ai, x) for ai in a]

V(a::Real, x::AbstractArray{<:Real}) = [V(a, xi) for xi in x]
V(a::AbstractArray{<:Real}, x::Real) = [V(ai, x) for ai in a]

dV(a::Real, x::AbstractArray{<:Real}) = [dV(a, xi) for xi in x]
dV(a::AbstractArray{<:Real}, x::Real) = [dV(ai, x) for ai in a]

W(a::Real, x::AbstractArray{<:Real}) = [W(a, xi) for xi in x]
W(a::AbstractArray{<:Real}, x::Real) = [W(ai, x) for ai in a]

dW(a::Real, x::AbstractArray{<:Real}) = [dW(a, xi) for xi in x]
dW(a::AbstractArray{<:Real}, x::Real) = [dW(ai, x) for ai in a]
