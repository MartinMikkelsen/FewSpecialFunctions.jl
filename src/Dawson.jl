# Adapted from the algorithms in:
# M. R. Zaghloul, Numerical Algorithms 95, 1291-1308 (2023).
# https://doi.org/10.1007/s11075-023-01608-8
# Source: https://github.com/mofrehzaghloul/Dawson (MIT)

const _DAWSON_MAXITER = 10_000

_dawson_tolerance(::Type{T}) where {T <: AbstractFloat} = T(4) * eps(T)
_dawson_small_start(::Type{T}) where {T <: AbstractFloat} = sqrt(eps(T))
_dawson_large_start(::Type{T}) where {T <: AbstractFloat} = max(T(8), sqrt(T(2) * log(inv(eps(T)))))
_dawson_asymptotic_start(::Type{T}) where {T <: AbstractFloat} = inv(sqrt(eps(T)))

function _dawson_small_cf_order(m::Int, x::T) where {T <: AbstractFloat}
    x2 = x * x
    sign = isodd(m) ? one(T) : -one(T)
    y = sign * T(2 * m) * x2 / T(4 * m^2 - 1)
    for k in (m - 1):-1:1
        sign = isodd(k) ? one(T) : -one(T)
        y = sign * T(2 * k) * x2 / (T(4 * k^2 - 1) * (one(T) + y))
    end
    return x / (one(T) + y)
end

function _dawson_small_cf(x::T) where {T <: AbstractFloat}
    previous = _dawson_small_cf_order(1, x)
    tol = _dawson_tolerance(T)
    for m in 2:_DAWSON_MAXITER
        current = _dawson_small_cf_order(m, x)
        if abs(current - previous) <= tol * max(one(T), abs(current))
            return current
        end
        previous = current
    end
    throw(ErrorException("Dawson small continued fraction failed to converge for x = $x"))
end

function _dawson_middle_series(x::T) where {T <: AbstractFloat}
    x2 = x * x
    sum = one(T)
    term = one(T)
    tol = _dawson_tolerance(T)
    for n in 0:(_DAWSON_MAXITER - 1)
        term *= x2 * T(2 * n + 1) / T((2 * n + 3) * (n + 1))
        sum += term
        if abs(term) <= tol * abs(sum)
            return x * exp(-x2) * sum
        end
    end
    throw(ErrorException("Dawson middle series failed to converge for x = $x"))
end

function _dawson_large_cf_order(m::Int, x::T) where {T <: AbstractFloat}
    half = inv(T(2))
    y = T(m) / x
    for k in (m - 1):-1:1
        y = T(k) / muladd(-half, y, x)
    end
    return inv(muladd(T(2), x, -y))
end

function _dawson_large_cf(x::T) where {T <: AbstractFloat}
    previous = _dawson_large_cf_order(1, x)
    tol = _dawson_tolerance(T)
    for m in 2:_DAWSON_MAXITER
        current = _dawson_large_cf_order(m, x)
        scale = abs(current)
        if scale == zero(T)
            if abs(current - previous) <= tol
                return current
            end
        elseif abs(current - previous) <= tol * scale
            return current
        end
        previous = current
    end
    throw(ErrorException("Dawson large continued fraction failed to converge for x = $x"))
end

function _dawson_real(x::T) where {T <: AbstractFloat}
    isnan(x) && return x
    iszero(x) && return x
    isinf(x) && return inv(x + x)
    ax = abs(x)
    value = if ax <= _dawson_small_start(T)
        _dawson_small_cf(ax)
    elseif ax <= _dawson_large_start(T)
        _dawson_middle_series(ax)
    elseif ax >= _dawson_asymptotic_start(T)
        inv(ax + ax)
    else
        _dawson_large_cf(ax)
    end
    return copysign(value, x)
end

dawson(x::T) where {T <: AbstractFloat} = _dawson_real(x)
dawson(x::Integer) = dawson(Float64(x))
dawson(x::Real) = dawson(float(x))
