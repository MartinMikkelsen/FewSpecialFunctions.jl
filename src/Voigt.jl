# Adapted from the Fourier-expansion method in:
# S. M. Abrarov, B. M. Quine, and R. K. Jagpal, “Rapidly convergent series for high-accuracy calculation of the Voigt function,” Journal of Quantitative Spectroscopy & Radiative Transfer 111, 372-375 (2010).
# https://doi.org/10.1016/j.jqsrt.2009.09.005

const _VOIGT_T_M = 12
const _VOIGT_FOURIER_COEFFICIENTS = (
    0.2954089751509193,
    0.27584023329217705,
    0.22457395522461585,
    0.1594149382739117,
    0.0986657664154542,
    0.053244140787639394,
    0.02505215000539365,
    0.01027746567053954,
    0.00367616433284485,
    0.001146493641242233,
    0.0003117570150461973,
    7.391433429603023e-5,
    1.527949342800837e-5,
    2.753956608221068e-6,
    4.327858781901254e-7,
    5.9300304087459056e-8,
    7.0844903077482305e-9,
    7.379520635816759e-10,
    6.70217160600201e-11,
    5.307265163470824e-12,
    3.664324113467629e-13,
    2.205894944941035e-14,
    1.1578268626285645e-15,
    5.2987114294673457e-17,
)

@inline function _voigt_integral(y::T, q::T, cos_qtM::T, sin_qtM::T) where {T <: AbstractFloat}
    tM = T(_VOIGT_T_M)
    # The direct form subtracts nearly equal terms when q*tM is small. The
    # expm1 form is stable throughout |q|*tM ≤ sqrt(eps(T)), including q = 0.
    abs(q) <= sqrt(eps(T)) / tM && return real(-expm1((-y + im * q) * tM) / (y - im * q))

    scale = max(abs(y), abs(q))
    y_scaled, q_scaled = y / scale, q / scale
    numerator = y_scaled - exp(-y * tM) * (y_scaled * cos_qtM - q_scaled * sin_qtM)
    return numerator / (scale * (y_scaled^2 + q_scaled^2))
end

function _voigt_fourier(x::T, y::T) where {T <: AbstractFloat}
    (isnan(x) || isnan(y)) && return T(NaN)
    (isinf(x) || isinf(y)) && return zero(T)

    tM = T(_VOIGT_T_M)
    sin_xtM, cos_xtM = sincos(x * tM)
    K = T(_VOIGT_FOURIER_COEFFICIENTS[1]) * _voigt_integral(y, x, cos_xtM, sin_xtM) / 2

    for n in eachindex(_VOIGT_FOURIER_COEFFICIENTS[2:end])
        nπ = T(n) * T(π) / tM
        sin_nπ, cos_nπ = sinpi(T(n)), cospi(T(n))
        sin_plus = sin_xtM * cos_nπ + cos_xtM * sin_nπ
        cos_plus = cos_xtM * cos_nπ - sin_xtM * sin_nπ
        sin_minus = sin_xtM * cos_nπ - cos_xtM * sin_nπ
        cos_minus = cos_xtM * cos_nπ + sin_xtM * sin_nπ
        K += T(_VOIGT_FOURIER_COEFFICIENTS[n + 1]) * (
            _voigt_integral(y, x + nπ, cos_plus, sin_plus) +
            _voigt_integral(y, x - nπ, cos_minus, sin_minus)
        ) / 2
    end

    return K / sqrt(T(π))
end

function _voigt_real(x::T, y::T) where {T <: AbstractFloat}
    y < zero(T) && throw(DomainError(y, "Voigt function requires y ≥ 0"))
    iszero(y) && return exp(-x^2)
    return _voigt_fourier(x, y)
end

voigt(x::T, y::T) where {T <: AbstractFloat} = _voigt_real(x, y)

"""
    voigt(x::Real, y::Real)

Evaluate the real Voigt function
```math
K(x, y) = \\frac{y}{\\pi} \\int_{-\\infty}^{\\infty}
    \\frac{e^{-t^2}}{(x - t)^2 + y^2} \\, \\mathrm{d}t,
\\qquad y \\ge 0.
```
For `y == 0`, `K(x, 0) = exp(-x^2)`. The implementation follows the
[Fourier-expansion method of Abrarov, Quine, and Jagpal](https://doi.org/10.1016/j.jqsrt.2009.09.005).
"""
function voigt(x::Real, y::Real)
    T = float(promote_type(typeof(x), typeof(y)))
    return voigt(T(x), T(y))
end
