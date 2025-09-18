using SpecialFunctions
export Clausen, Ci_complex, f_n, F_clausen

const ξ10 = [
    1.02058363572825669e-1,
    1.10087335344093285e0,
    3.87538426676163313e0,
    9.49089900532915982e0,
    1.92369217680950509e1,
    3.48489807870071239e1,
    5.88422971202513756e1,
    9.52057687864990326e1,
    1.51245498588511132e2,
    2.44769266678480451e2,
]

const A10 = [
    1.03011687962788351e0,
    2.2580220852505727e-1,
    1.67040618192678355e-2,
    6.0460703772163201e-4,
    1.16718411051188568e-5,
    1.1536011315284419e-7,
    5.23125205372251291e-10,
    8.88585842372452499e-13,
    3.78823474931504689e-16,
    1.45749170427449731e-20,
]
const ξ20 = [
    7.8332741180217739859429592034215121e-2,
    7.99480557580280868888504305965764455e-1,
    2.633559843707676053063024600681851e0,
    6.067125354800522713778706639849939e0,
    1.1578014308014972867250241260973853e1,
    1.967574219682007399786022325652155e1,
    3.09284092373144262920338413642651e1,
    4.5988905145672958078531429470368746e1,
    6.5602547968040040431325251181403039e1,
    9.0652159037294131940072043099536872e1,
    1.2218773346692348159267034032123833e2,
    1.6148763659034474884235758941294436e2,
    2.1014062487874825104404326118075493e2,
    2.7017538454329128248632034154940478e2,
    3.4427146853197803740611415515669254e2,
    4.3612864159248546858458013080006752e2,
    5.5117670639883226012194621182330183e2,
    6.9813208096151573112006881397004676e2,
    8.9318894895282773667447559458923804e2,
    1.1765731082977506193890547798997328e3,
]

const A20 = [
    9.4521844770213813522739498751084446e-1,
    2.8673730910005180272566734208681658e-1,
    3.7990715886134892035536939814054822e-2,
    3.1053688857936936821748901103472756e-3,
    1.7979523985634735267965413590172828e-4,
    7.6599061254187375464405113636273785e-6,
    2.4224398021113042502544194520811499e-7,
    5.67246161661686603385118458797347e-9,
    9.740052786086906643974014704731242e-11,
    1.2068824524541185034179793734503826e-12,
    1.0549929160747159130056637444912904e-14,
    6.3107861649521708824465976542847217e-17,
    2.4806740961948186082790889707150879e-19,
    6.068944075643647849080969319008971e-22,
    8.5770997286711710319485672039762078e-25,
    6.296750414836763110694582719147875e-28,
    2.0455706071903958823179960004147795e-31,
    2.258006419692752335532584347360038e-35,
    5.1214712907227134883116239441893341e-40,
    6.6554055933455521526556015293110494e-46,
]
"""
    f_n(n::Int, k::Int, θ::Float64)

Compute the Clausen series summand fₙ(k, θ):
    sin(kθ)/kⁿ for even n,
    cos(kθ)/kⁿ for odd n.
"""
function f_n(n::Int, k::Int, θ::Real)
    θ_float = Float64(θ)

    if iseven(n)
        return sin(k * θ_float) / (float(k))^n
    else
        return cos(k * θ_float) / (float(k))^n
    end
end

"""
    F_clausen(n::Int, z::ComplexF64, θ::Float64)

Dispatch to the correct primitive function Fₙ(z, θ) for n = 1..6.
"""
function F_clausen(n::Int, z::ComplexF64, θ::Float64)
    if n == 1
        return F1(z, θ)
    elseif n == 2
        return F2(z, θ)
    elseif n == 3
        return F3(z, θ)
    elseif n == 4
        return F4(z, θ)
    elseif n == 5
        return F5(z, θ)
    elseif n == 6
        return F6(z, θ)
    else
        throw(ArgumentError("Only Fₙ with n = 1 to 6 are implemented."))
    end
end

"""
    Ci_complex(z::ComplexF64)

Complex cosine integral function used in Clausen function calculations.
"""
function Ci_complex(z::ComplexF64)
    if isinf(real(z)) && imag(z) == 0
        return real(z) > 0 ? (0.0 + 0im) : (π * im)
    end

    if isinf(z)
        return NaN + NaN * im
    end

    if z == 0.0 + 0.0im
        return NaN + NaN * im
    end

    v = -0.5 * (expint(im * z) + expint(-im * z))

    if real(z) < 0
        v += π * im
    end

    return ComplexF64(v)
end

const γ = Base.MathConstants.γ

@inline F1(z::ComplexF64, θ::Float64) = Ci_complex(z * θ)

@inline F2(z::ComplexF64, θ::Float64) =
    (θ * z * Ci_complex(z * θ) - sin(θ * z)) / z

@inline F3(z::ComplexF64, θ::Float64) =
    -1 / (2 * z^2) * (
    θ^2 * z^2 * Ci_complex(z * θ) + cos(θ * z) - θ * z * sin(θ * z)
)

@inline F4(z::ComplexF64, θ::Float64) =
    -1 / (6 * z^3) * (
    θ^3 * z^3 * Ci_complex(z * θ) + (2 - θ^2 * z^2) * sin(θ * z) + θ * z * cos(θ * z)
)

@inline F5(z::ComplexF64, θ::Float64) =
    1 / (24 * z^4) * (
    θ^4 * z^4 * Ci_complex(z * θ)
        + θ * z * (2 - θ^2 * z^2) * sin(θ * z)
        + (θ^2 * z^2 - 6) * cos(θ * z)
)


@inline F6(z::ComplexF64, θ::Float64) =
    1 / (120 * z^5) * (
    θ^5 * z^5 * Ci_complex(z * θ)
        + θ * z * (θ^2 * z^2 - 6) * cos(θ * z) - (θ^4 * z^4 - 2 * θ^2 * z^2 + 24) * sin(θ * z)
)

"""
    Clausen(n::Int, θ::Float64; N::Int=10, m::Int=20)

Compute the Clausen function of order `n` at angle `θ`.

References:
- [Clausen function](https://en.wikipedia.org/wiki/Clausen_function)
- [Implementation paper](https://doi.org/10.1007/s10543-023-00944-4)
"""
function Clausen(n::Int, θ::Float64; N::Int = 10, m::Int = 20)
    n < 1 || n > 6 && throw(ArgumentError("Only n=1..6 supported"))
    (N == 10 || N == 20) || throw(ArgumentError("Only N=10 or N=20 supported"))

    if θ == 0.0
        return iseven(n) ? 0.0 : zeta(n)
    end

    θmod = mod(θ, 2π)

    if θmod <= π
        φ, sign = θmod, 1.0
    else
        φ, sign = 2π - θmod, (-1.0)^(n + 1)
    end

    if n == 1
        if isapprox(φ, 0.0; atol = 1.0e-14) || isapprox(φ, 2π; atol = 1.0e-14)
            return Inf
        end
        return sign * (-log(abs(2 * sin(φ / 2))))
    end

    ξ, A = N == 10 ? (ξ10, A10) : (ξ20, A20)

    S1 = sum(f_n(n, k, φ) for k in 1:(m - 1))

    S2 = 0.0
    for ν in 1:N
        z = (m - 0.5) + 0.5im * sqrt(ξ[ν])
        S2 += A[ν] * real(F_clausen(n, z, φ))
    end

    result0 = S1 - (π / 4) * S2
    return sign * result0
end
