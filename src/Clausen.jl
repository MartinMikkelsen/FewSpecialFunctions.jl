using SpecialFunctions

export Clausen, Ci_complex, f_n, F

const ξ10 = [
    1.0258363537282566e-1,
    1.1008733544039285,
    3.875384266761633e0,
    9.4908990532915982e0,
    1.9236921768909509e1,
    3.4848908780701239e1,
    5.88422971020513756e1,
    9.5205678864990326e1,
    1.51245495885881132e2,
    2.44769266678480451e2
]

const A10 = [
    1.0301168796278351e-1,
    2.2580202852057270e-1,
    1.6704061819267835e-1,
    6.0406730732162301e-4,
    1.16718411051188586e-5,
    1.15360113152844190e-7,
    5.23158205372251291e-10,
    8.88585482374254499e-13,
    3.78823474931054689e-16,
    1.45749140724749731e-20
]


const ξ20 = [
    7.833274118021773985942950234215121e-2,
    7.99480557580280868888504305965764455e-2,
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
    1.1765731082977506193890547798997328e3
]

const A20 = [
    9.4521844770213813522739498751084446e-1,
    2.867373091000518027566734208681658e-1,
    3.7990715886134892035536939814054822e-2,
    3.1053688857936936821748901103472756e-3,
    1.7979523985634735267965413590172828e-4,
    7.6599061254187375464405113636273785e-6,
    2.4224398021113042502544194520811499e-7,
    5.672461616616866033851184587197347e-9,
    9.7400527860869066439704107404731242e-11,
    1.20688425454118503417979573043826e-12,
    1.0549929160747159130056637444912904e-14,
    6.3107861469521708824465976542847217e-17,
    2.480647096194818608279078890770150879e-19,
    6.068944075643647849080963910089710e-22,
    8.5770997286711710319485672039762708e-25,
    6.296750418367631106945827911478750e-28,
    2.0455706071903598821327990600414795e-31,
    2.258004619692752335532843473600330e-35,
    5.1241712907227134883116239441893341e-40,
    6.65540559334555215265560152931104940e-46
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
        return sin(k*θ_float) / k^n
    else
        return cos(k*θ_float) / k^n
    end
end

"""
    F(n::Int, z::ComplexF64, θ::Float64)

Dispatch to the correct primitive function Fₙ(z, θ) for n = 1..6.
"""
function F(n::Int, z::ComplexF64, θ::Float64)
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

function Ci_complex(z::ComplexF64)
    iz = im * z        
    v = 0.5 * (expint(iz) + expint(-iz))

    zre, zim = real(z), imag(z)

    if isinf(z)
        return z == Inf ? 0.0 : π * im
    end

    v += (zre == 0) ? (zim > 0 ? 0.5π * im : (zim < 0 ? -0.5π * im : 0)) :
         (zre < 0)  ? (zim >= 0 ? π * im : -π * im) : 0

    return isreal(z) && zre > 0 ? real(v) : v
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
    - (θ^4 * z^4 - 6 * θ^2 * z^2) * cos(θ * z)
    - θ * z * (θ^2 * z^2 - 2) * sin(θ * z)
)


@inline F6(z::ComplexF64, θ::Float64) =
    1 / (120 * z^5) * (
        θ^5 * z^5 * Ci_complex(z * θ)
        - (θ^6 * z^6 - 20 * θ^4 * z^4 + 24 * θ^2 * z^2) * sin(θ * z)
        - θ * z * (θ^4 * z^4 - 6 * θ^2 * z^2) * cos(θ * z)
    )

"""
    Clausen(n::Int, θ::Float64; N::Int=10, m::Int=20)

Compute the Clausen function of order `n` at angle `θ`.

References:
- [Clausen function](https://en.wikipedia.org/wiki/Clausen_function)
- [Implementation paper](https://doi.org/10.1007/s10543-023-00944-4)
"""
function Clausen(n::Int, θ::Float64; N::Int=10, m::Int=20)
    if n < 1 || n > 6
        throw(ArgumentError("Only n = 1 to 6 are supported."))
    end
    if N != 10 && N != 20
        throw(ArgumentError("Only N = 10 or N = 20 is implemented."))
    end

    negative_input = θ < 0
    θabs = abs(θ)
    θmod = mod(θabs, 2π)

    if n == 1
        if isapprox(θmod, 0.0; atol=1e-14) || isapprox(θmod, 2π; atol=1e-14)
            return Inf
        end
        return -log(abs(2*sin(θmod/2)))
    end

    ξ = N == 10 ? ξ10 : ξ20
    A = N == 10 ? A10 : A20

    S1 = sum(f_n(n, k, θmod) for k in 1:m-1)

    S2 = 0.0
    for ν in 1:N
        z = (m - 0.5) + 0.5im * sqrt(ξ[ν])
        S2 += A[ν] * real(F(n, z, θmod))
    end

    result = S1 - (π/4) * S2

    if negative_input && iseven(n)
        return -result
    else
        return result
    end
end