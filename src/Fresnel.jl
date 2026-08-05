# Adapted from the algorithms in:
# M. R. Zaghloul and L. Alrawas, Numerical Algorithms 96, 489-506 (2024).
# https://doi.org/10.1007/s11075-023-01654-2
# Source: https://github.com/mofrehzaghloul/Fresnel_Integrals
# Copyright (c) 2023 Mofreh R Zaghloul and Leen Alrawas
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

using SpecialFunctions

export fresnel, FresnelS, FresnelC, FresnelE

function _fresnel_series(z)
    T = typeof(real(z))
    u = (T(π) / 2) * z^2
    u2 = u^2
    C = z
    S = z * u / 3
    Cterm = C
    Sterm = S
    tolerance = 4eps(T)

    for k in 1:10_000
        Cterm *= -((4k - 3) * u2) / (2k * (2k - 1) * (4k + 1))
        Sterm *= -((4k - 1) * u2) / (2k * (2k + 1) * (4k + 3))
        C += Cterm
        S += Sterm

        max(abs(Cterm), abs(Sterm)) <= tolerance * max(one(T), abs(C), abs(S)) && break
    end

    return C, S
end

function _fresnel_miller(z, parity, fnext)
    T = typeof(real(z))
    u = (T(π) / 2) * z^2
    K = max(32, ceil(Int, 2abs(u) + precision(T)))
    f = one(z)
    total = zero(z)

    for k in K:-1:0
        fk = (2k + 3) * f / u - fnext
        k % 2 == parity && (total += fk)
        fnext, f = f, fk
    end

    return total * sin(u) * z / (f * u)
end

function _fresnel_asymptotic_auxiliary(z, offset)
    T = typeof(real(z))
    πT = T(π)
    z4 = z^4
    term = one(z)
    sum = term
    tolerance = 4eps(T)

    for k in 1:10_000
        nextterm = -term * (4k - 3 + offset) * (4k - 1 + offset) / (πT^2 * z4)
        abs(nextterm) > abs(term) && break
        sum += nextterm
        abs(nextterm) <= tolerance * abs(sum) && break
        term = nextterm
    end

    return sum
end

function _fresnel_asymptotic(z)
    T = typeof(real(z))
    πT = T(π)
    θ = (πT / 2) * z^2
    f = _fresnel_asymptotic_auxiliary(z, 0) / (πT * z)
    g = _fresnel_asymptotic_auxiliary(z, 2) / (πT^2 * z^3)

    if abs(imag(z)) > abs(real(z))
        exponential = exp(im * θ) / 2
        f += (one(z) + im) * exponential
        g += (one(z) - im) * exponential
    end

    S = one(z) / 2 - f * cos(θ) - g * sin(θ)
    C = one(z) / 2 + f * sin(θ) - g * cos(θ)
    return C, S
end

function _fresnel_erf_real(x::T) where {T <: AbstractFloat}
    w = (sqrt(T(π)) / 2) * (one(T) - im) * x
    E = (one(T) + im) / 2 * erf(w)
    return real(E), imag(E)
end

function _fresnel_asymptotic_start(::Type{T}) where {T <: AbstractFloat}
    return max(T(5), sqrt(2 * log(inv(eps(T))) / T(π)) + one(T))
end

function _fresnel_pair(z)
    T = typeof(real(z))
    z2 = abs2(z)

    if z2 <= T(6.9)
        return _fresnel_series(z)
    elseif z2 <= T(25 + 2precision(T))
        C = _fresnel_miller(z, 0, zero(z))
        S = _fresnel_miller(z, 1, one(z))
        return C, S
    else
        return _fresnel_asymptotic(z)
    end
end

function _fresnel_real(x::T) where {T <: AbstractFloat}
    if iszero(x)
        return zero(T), zero(T), zero(Complex{T})
    end

    ax = abs(x)
    C, S = if ax^2 <= T(6.9)
        _fresnel_series(ax)
    elseif ax < _fresnel_asymptotic_start(T)
        _fresnel_erf_real(ax)
    else
        _fresnel_asymptotic(ax)
    end

    if x < zero(T)
        C = -C
        S = -S
    end

    return C, S, complex(C, S)
end

function _fresnel_complex(z::Complex{T}) where {T <: AbstractFloat}
    x = real(z)
    y = imag(z)

    if iszero(y)
        C, S, E = _fresnel_real(x)
        return complex(C), complex(S), E
    elseif iszero(x)
        C, S, _ = _fresnel_real(y)
        return im * C, -im * S, im * C + S
    elseif x < zero(T)
        C, S, E = _fresnel_complex(-z)
        return -C, -S, -E
    end

    C, S = _fresnel_pair(z)
    return C, S, C + im * S
end

"""
    fresnel(z::Number)

Compute the Fresnel cosine and sine integrals and their combination:

    S(z) = ∫₀ᶻ sin(π / 2 * t²) dt
    C(z) = ∫₀ᶻ cos(π / 2 * t²) dt

Returns `(C, S, C + im * S)`. For complex inputs, `C` and `S` are the analytic
complex Fresnel integrals.
"""
fresnel(z::T) where {T <: AbstractFloat} = _fresnel_real(z)
fresnel(z::Integer) = fresnel(Float64(z))
fresnel(z::Complex{T}) where {T <: AbstractFloat} = _fresnel_complex(z)
fresnel(z::Complex{<:Integer}) = fresnel(ComplexF64(z))
fresnel(z::Real) = fresnel(float(z))
fresnel(z::Complex) = fresnel(complex(float(real(z)), float(imag(z))))

"""
    FresnelC(z::Number)

Compute the Fresnel cosine integral `C(z)`.
"""
FresnelC(z::Number) = fresnel(z)[1]

"""
    FresnelS(z::Number)

Compute the Fresnel sine integral `S(z)`.
"""
FresnelS(z::Number) = fresnel(z)[2]

"""
    FresnelE(z::Number)

Compute the Fresnel integral combination `C(z) + im * S(z)`.
"""
FresnelE(z::Number) = fresnel(z)[3]
