# DLMF 13.14.2–3, 13.14.33, 13.19.3; see also Thompson & Barnett
# (1986), sections 2–3. Reuse the Kummer series, with extra precision for
# cancellation and large prefactors. No ODE solver dependency is needed.

_whittaker_gamma(z) = isreal(z) ? gamma(real(z)) : exp(loggamma(z))
_whittaker_rgamma(z) = isreal(z) && real(z) <= 0 && isinteger(real(z)) ? zero(z) : inv(_whittaker_gamma(z))

function _whittaker_m(κ, μ, z)
    return exp(-z / 2 + (μ + 1 / 2) * log(z)) * HypergeometricFunctions._₁F₁(μ - κ + 1 / 2, 2μ + 1, z)
end

function _whittaker_w_connection(κ, μ, z)
    return _whittaker_gamma(-2μ) * _whittaker_rgamma(1 / 2 - μ - κ) * _whittaker_m(κ, μ, z) +
        _whittaker_gamma(2μ) * _whittaker_rgamma(1 / 2 + μ - κ) * _whittaker_m(κ, -μ, z)
end

function _whittaker_w_series(κ, μ, z, p)
    if isreal(μ) && isinteger(2real(μ))
        # The connection formula has removable poles at integer 2μ. A
        # symmetric limit has O(h²) error; reserve bits for its cancellation.
        h = ldexp(one(BigFloat), -p - 16)
        return (_whittaker_w_connection(κ, μ + h, z) + _whittaker_w_connection(κ, μ - h, z)) / 2
    end
    return _whittaker_w_connection(κ, μ, z)
end

# Only accept the Poincaré expansion when decreasing terms reach precision.
# This also evaluates terminating polynomial cases without connection poles.
function _whittaker_w_asymptotic(κ, μ, z)
    r = s = one(z)
    ds = zero(z)
    for n in 1:1000
        next = -r * ((μ - κ + 1 / 2 + (n - 1)) / z) * (-μ - κ + 1 / 2 + (n - 1)) / n
        abs(next) > abs(r) && return nothing
        s += next
        ds -= n * next / z
        r = next
        if abs(r) <= eps(real(z)) * abs(s) / 16
            scale = exp(-z / 2 + κ * log(z))
            return scale * s, scale * (ds + (κ / z - 1 / 2) * s)
        end
    end
    return nothing
end

function _whittaker(κ::Number, μ::Number, z::Number, kind::Symbol, derivative::Bool)
    κ, μ, z = promote(float(κ), float(μ), float(z))
    all(isfinite, (κ, μ, z)) || throw(DomainError((κ, μ, z), "Whittaker parameters must be finite"))
    iszero(z) && throw(DomainError(z, "Whittaker functions require a nonzero argument"))
    z isa Real && z < 0 && throw(DomainError(z, "use a complex argument to select a Whittaker branch at negative z"))
    kind === :M && isreal(μ) && 2real(μ) <= -1 && isinteger(2real(μ)) &&
        throw(DomainError(μ, "Whittaker M has a parameter pole at negative integer 2μ"))
    T = typeof(z)
    p = precision(real(z))
    # Share the MPFR lock with cylinder fallbacks on Julia 1.10.
    result = lock(_cylinder_precision_lock) do
        if kind === :W
            asymptotic = setprecision(BigFloat, p + 32) do
                _whittaker_w_asymptotic(big(κ), big(μ), big(z))
            end
            if asymptotic !== nothing
                return derivative ? last(asymptotic) : first(asymptotic)
            end
        end
        # Two exponential scales can cancel in W's connection formula.
        # Reserve bits for growth and proximity to its removable poles.
        extra = ceil(Int, (2abs(z) + 2sqrt(abs((μ - κ + 1 / 2) * z)) + abs(κ) + abs(μ)) / log(2)) + 64
        distance = abs(2μ - round(2real(μ)))
        iszero(distance) || (extra += max(0, -exponent(distance)))
        setprecision(BigFloat, 2p + extra) do
            kb, mb, zb = big(κ), big(μ), big(z)
            y = kind === :M ? _whittaker_m(kb, mb, zb) : _whittaker_w_series(kb, mb, zb, p)
            derivative || return y
            a = mb - kb + 1 / 2
            shifted = kind === :M ? a / (2mb + 1) * _whittaker_m(kb - 1 / 2, mb + 1 / 2, zb) :
                -a * _whittaker_w_series(kb - 1 / 2, mb + 1 / 2, zb, p)
            return (-1 / 2 + (mb + 1 / 2) / zb) * y + shifted / sqrt(zb)
        end
    end
    if T <: Real
        return T === BigFloat ? BigFloat(real(result); precision = p) : T(real(result))
    end
    return T === Complex{BigFloat} ? complex(BigFloat(real(result); precision = p), BigFloat(imag(result); precision = p)) : T(result)
end

"""
    WhittakerM(κ::Number, μ::Number, z::Number)

Compute the Whittaker function `Mκμ(z) = exp(-z/2) * z^(μ+1/2) *
₁F₁(μ-κ+1/2, 1+2μ, z)` on the principal branch.

Inputs must be finite and `z ≠ 0`; negative real arguments require an
explicit complex input. Negative integer `2μ` is a parameter pole.
Invalid inputs raise `DomainError`. Real inputs with `z > 0` return real
values; complex inputs return complex values. Float32, Float64 and BigFloat
are supported, with mixed inputs promoted. Use broadcasting for arrays.

Evaluation uses a precision-guarded Kummer series. Extra precision controls
cancellation and intermediate overflow, at a cost in speed for large inputs.
See [DLMF 13.14](https://dlmf.nist.gov/13.14) and Thompson & Barnett (1986),
[doi:10.1016/0021-9991(86)90046-X](https://doi.org/10.1016/0021-9991(86)90046-X).
"""
WhittakerM(κ::Number, μ::Number, z::Number) = _whittaker(κ, μ, z, :M, false)

"""
    WhittakerW(κ::Number, μ::Number, z::Number)

Compute `Wκμ(z) = exp(-z/2) * z^(μ+1/2) * U(μ-κ+1/2, 1+2μ, z)`,
where `U` here is Tricomi's confluent hypergeometric function.

The domain, principal-branch convention and type behavior match
[`WhittakerM`](@ref), but `W` also accepts negative integer `2μ`.
Uses an asymptotic expansion when it reaches working precision, otherwise
a guarded connection formula, including limits at integer `2μ`.
These are series methods described by Thompson & Barnett (1986); this is
not a port of their full COULCC algorithm.
"""
WhittakerW(κ::Number, μ::Number, z::Number) = _whittaker(κ, μ, z, :W, false)

"""
    dWhittakerM(κ::Number, μ::Number, z::Number)

Compute the argument derivative of [`WhittakerM`](@ref) using the analytic
derivative of its Kummer representation. The domain and type behavior match
`WhittakerM`.
"""
dWhittakerM(κ::Number, μ::Number, z::Number) = _whittaker(κ, μ, z, :M, true)

"""
    dWhittakerW(κ::Number, μ::Number, z::Number)

Compute the argument derivative of [`WhittakerW`](@ref) using its asymptotic
expansion or the analytic derivative of its Tricomi representation. The
domain and type behavior match `WhittakerW`.
"""
dWhittakerW(κ::Number, μ::Number, z::Number) = _whittaker(κ, μ, z, :W, true)
