using SpecialFunctions

export FermiDiracIntegral, FermiDiracIntegralNorm

@doc raw"""
    FermiDiracIntegral(j::Real, x::Real)

Approximate the unnormalized complete Fermi–Dirac integral

```math
\mathcal{F}_j(x) = \int_0^\infty \frac{t^j}{\exp(t-x)+1}\,\mathrm{d}t
```

of order `j ≥ -1/2` at real `x`. No factor `1/Γ(j+1)` is included; see
[`FermiDiracIntegralNorm`](@ref) for the normalized form. Orders `j < -1/2`
throw an `ErrorException`.

The arguments are promoted to a common floating-point type `T`, which is also
the return type. The method then depends on the order:

| Order `j` | Method | Measured relative error (`Float64`) |
|:--|:--|:--|
| `0` | closed form `log(1 + exp(x))` | `1.5e-14` for `x ≥ -5`; see below |
| `-1/2`, `1/2`, `3/2`, `5/2` | rational approximations of Antia [5]: in `exp(x)` for `x < 2`, in `1/x²` for `x ≥ 2` | `2.8e-12`, `5.4e-13`, `5.1e-13`, `2.5e-13` |
| any other `j > -1/2` | closed-form expression of Aymerich-Humet et al. [4] | `5.8e-3` to `1.2e-2` at the tested orders `1/4`, `1`, `2`, `3`, `9/2`; `2.9e-2` at order `10` |

Each figure is the maximum over a grid of 190 points in `-30 ≤ x ≤ 60`,
relative to quadrature of the integral in 512-bit arithmetic. They are
measurements at the listed orders, not bounds; `docs/accuracy_checks.jl`
reproduces them. The order is matched by exact comparison after conversion to
`T`, so `1/2` selects the rational approximation and `0.5 + 1e-12` does not.

The return type does not change the method. The coefficients of [5] are
double-precision constants and the expression of [4] is an approximation, so a
`BigFloat` result has the same error as the `Float64` one, and a `Float32`
result is limited by single precision. Two limits of the `j = 0` form: for
`x ≪ 0` the sum `1 + exp(x)` rounds, so the relative error grows as `x`
decreases (about `1e-3` at `x = -30` in `Float64`) although the absolute
error stays below `eps(T)`, and for `x` large enough that `exp(x)` overflows
(`x > 709` in `Float64`) the result is `Inf`.

# Examples
```jldoctest
julia> FermiDiracIntegral(3 / 2, 1.0) ≈ 2.6616826247307124
true

julia> FermiDiracIntegral(0, 0.0) ≈ log(2)
true

julia> FermiDiracIntegral(0.5f0, 1.0f0) isa Float32
true
```

# References
1. D. Bednarczyk and J. Bednarczyk, Phys. Lett. A 64, 409 (1978)
2. J. S. Blakemore, Solid-State Electron. 25, 1067 (1982)
3. X. Aymerich-Humet, F. Serra-Mestres, and J. Millan, Solid-State Electron. 24, 981 (1981)
4. X. Aymerich-Humet, F. Serra-Mestres, and J. Millan, J. Appl. Phys. 54, 2850 (1983)
5. H. M. Antia, Astrophys. J. Suppl. Ser. 84, 101 (1993)

See also [Kim and Lundstrom, *Notes on Fermi-Dirac Integrals*](https://arxiv.org/abs/0811.0116)
and [DLMF 25.12(iii)](https://dlmf.nist.gov/25.12#iii).
"""
function FermiDiracIntegral(j::Real, x::Real)
    T = float(promote_type(typeof(j), typeof(x)))
    jT = T(j)
    xT = T(x)
    return _FermiDiracIntegral(jT, xT)
end

function _FermiDiracIntegral(j::T, x::T) where {T <: AbstractFloat}
    if j < -one(T) / 2
        error("The order should be equal to or larger than -1/2.")
    elseif j == 0
        y = log(one(T) + exp(x))
    elseif j == -one(T) / 2
        # Method from [5]
        a1 = (T(1.71446374704454e7), T(3.88148302324068e7), T(3.16743385304962e7), T(1.14587609192151e7), T(1.83696370756153e6), T(1.14980998186874e5), T(1.98276889924768e3), one(T))
        b1 = (T(9.67282587452899e6), T(2.87386436731785e7), T(3.26070130734158e7), T(1.77657027846367e7), T(4.81648022267831e6), T(6.13709569333207e5), T(3.13595854332114e4), T(4.35061725080755e2))
        a2 = (T(-4.46620341924942e-15), T(-1.58654991146236e-12), T(-4.44467627042232e-10), T(-6.84738791621745e-8), T(-6.64932238528105e-6), T(-3.69976170193942e-4), T(-1.12295393687006e-2), T(-1.60926102124442e-1), T(-8.52408612877447e-1), T(-7.45519953763928e-1), T(2.98435207466372e0), one(T))
        b2 = (T(-2.23310170962369e-15), T(-7.94193282071464e-13), T(-2.22564376956228e-10), T(-3.43299431079845e-8), T(-3.33919612678907e-6), T(-1.86432212187088e-4), T(-5.69764436880529e-3), T(-8.34904593067194e-2), T(-4.7877084400944e-1), T(-4.99759250374148e-1), T(1.86795964993052e0), T(4.16485970495288e-1))
        y = _fermi_dirac_rational(j, x, a1, b1, a2, b2)
    elseif j == one(T) / 2
        # Method from [5]
        a1 = (T(5.75834152995465e6), T(1.30964880355883e7), T(1.07608632249013e7), T(3.93536421893014e6), T(6.4249323371564e5), T(4.16031909245777e4), T(7.77238678539648e2), one(T))
        b1 = (T(6.49759261942269e6), T(1.70750501625775e7), T(1.6928813485616e7), T(7.95192647756086e6), T(1.83167424554505e6), T(1.95155948326832e5), T(8.17922106644547e3), T(9.02129136642157e1))
        a2 = (T(4.85378381173415e-14), T(1.64429113030738e-11), T(3.76794942277806e-9), T(4.69233883900644e-7), T(3.40679845803144e-5), T(1.32212995937796e-3), T(2.60768398973913e-2), T(2.48653216266227e-1), T(1.08037861921488e0), T(1.91247528779676e0), one(T))
        b2 = (T(7.28067571760518e-14), T(2.45745452167585e-11), T(5.62152894375277e-9), T(6.96888634549649e-7), T(5.02360015186394e-5), T(1.92040136756592e-3), T(3.66887808002874e-2), T(3.24095226486468e-1), T(1.16434871200131e0), T(1.34981244060549e0), T(2.0131183697593e-1), T(-2.14562434782759e-2))
        y = _fermi_dirac_rational(j, x, a1, b1, a2, b2)
    elseif j == 3 * one(T) / 2
        # Method from [5]
        a1 = (T(4.32326386604283e4), T(8.55472308218786e4), T(5.95275291210962e4), T(1.77294861572005e4), T(2.2187660779646e3), T(9.90562948053193e1), one(T))
        b1 = (T(3.25218725353467e4), T(7.01022511904373e4), T(5.50859144223638e4), T(1.959420745764e4), T(3.20803912586318e3), T(2.20853967067789e2), T(5.05580641737527e0), T(1.99507945223266e-2))
        a2 = (T(2.80452693148553e-13), T(8.60096863656367e-11), T(1.62974620742993e-8), T(1.6359884375205e-6), T(9.12915407846722e-5), T(2.62988766922117e-3), T(3.85682997219346e-2), T(2.78383256609605e-1), T(9.02250179334496e-1), one(T))
        b2 = (T(7.01131732871184e-13), T(2.10699282897576e-10), T(3.94452010378723e-8), T(3.84703231868724e-6), T(2.04569943213216e-4), T(5.31999109566385e-3), T(6.39899717779153e-2), T(3.14236143831882e-1), T(4.70252591891375e-1), T(-2.15540156936373e-2), T(2.34829436438087e-3))
        y = _fermi_dirac_rational(j, x, a1, b1, a2, b2)
    elseif j == 5 * one(T) / 2
        # Method from [5]
        a1 = (T(6.61606300631656e4), T(1.20132462801652e5), T(7.67255995316812e4), T(2.10427138842443e4), T(2.44325236813275e3), T(1.02589947781696e2), one(T))
        b1 = (T(1.99078071053871e4), T(3.79076097261066e4), T(2.60117136841197e4), T(7.97584657659364e3), T(1.10886130159658e3), T(6.35483623268093e1), T(1.16951072617142e0), T(3.31482978240026e-3))
        a2 = (T(8.42667076131315e-12), T(2.31618876821567e-9), T(3.54323824923987e-7), T(2.77981736000034e-5), T(1.14008027400645e-3), T(2.32779790773633e-2), T(2.39564845938301e-1), T(1.24415366126179e0), T(3.18831203950106e0), T(3.42040216997894e0), one(T))
        b2 = (T(2.94933476646033e-11), T(7.68215783076936e-9), T(1.12919616415947e-6), T(8.09451165406274e-5), T(2.81111224925648e-3), T(3.99937801931919e-2), T(2.27132567866839e-1), T(5.3188604522268e-1), T(3.70866321410385e-1), T(2.27326643192516e-2))
        y = _fermi_dirac_rational(j, x, a1, b1, a2, b2)
    else
        # Model proposed in [4]
        # Expressions from eqs. (6)-(7) of [4]
        a = (one(T) + T(15) / 4 * (j + one(T)) + (j + one(T))^2 / T(40))^(one(T) / 2)
        b = T(1.8) + T(0.61) * j
        c = 2 + (2 - sqrt(T(2))) * T(2)^(-j)
        y = (
            (j + one(T)) * T(2)^(j + one(T)) / (b + x + (abs(x - b)^c + a^c)^(one(T) / c))^(j + one(T)) +
                exp(-x) / gamma(j + one(T))
        )^-1
    end

    return y
end

function _fermi_dirac_rational(j::T, x::T, a1, b1, a2, b2) where {T <: AbstractFloat}
    if x < T(2)
        xx = exp(x)
        num = @evalpoly(xx, a1...)
        den = @evalpoly(xx, b1...)
        y = xx * num / den
    else
        xx = one(T) / x^2
        num = @evalpoly(xx, a2...)
        den = @evalpoly(xx, b2...)
        y = x^(j + one(T)) * num / den
    end

    return y
end

@doc raw"""
    FermiDiracIntegralNorm(j::Real, x::Real)

Approximate the normalized complete Fermi–Dirac integral

```math
F_j(x) = \frac{1}{\Gamma(j+1)}\int_0^\infty \frac{t^j}{\exp(t-x)+1}\,\mathrm{d}t
```

of order `j ≥ -1/2` at real `x`, computed as
`FermiDiracIntegral(j, x) / gamma(j + 1)`. With this normalization
`F_j(x) → exp(x)` as `x → -∞` and `dF_j/dx = F_{j-1}`.

The order restriction, the choice of method for each order, and the relative
error are those of [`FermiDiracIntegral`](@ref). In particular, only
`j ∈ (-1/2, 1/2, 3/2, 5/2)` and `j = 0` are computed to better than the
relative error of `6e-3` to `3e-2` measured for the other tested orders.

# Examples
```jldoctest
julia> FermiDiracIntegralNorm(1 / 2, 1.0) ≈ FermiDiracIntegral(1 / 2, 1.0) / (sqrt(π) / 2)
true

julia> FermiDiracIntegralNorm(0, 0.0) ≈ log(2)
true
```
"""
FermiDiracIntegralNorm(j::Real, eta::Real) = FermiDiracIntegral(j, eta) / gamma(j + 1)
