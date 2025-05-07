using QuanticsTCI
using SpecialFunctions
using Base: allequal

export MarcumQ

Q(μ, x) = gamma_inc(μ, x)[2]

# continued‐fraction ratio for scaled I_μ
c_μ(μ, ξ) = besseli(μ, ξ) / besseli(μ-1, ξ)

# log of Aₙ from eq. (32)
lnA(n, μ) = loggamma(μ + 0.5 + n) - loggamma(μ + 0.5 - n) - n*log(2) - loggamma(n+1)

# ζ²/2 from eq. (84)
function half_ζ2(x, y)
    δ = y - x - 1
    if abs(δ) < 1e-3
        c = [1.0,
             -(3*x + 1)/3,
             ((72*x + 42)*x + 7)/36,
             -(((2700*x + 2142)*x + 657)*x + 73)/540,
             ((((181440*x + 177552)*x + 76356)*x + 15972)*x + 1331)/12960]
        z = δ / (2*x + 1)^2
        s = sum(c[k] * z^k for k in 1:length(c))
        return (2*x + 1)^3 * s^2 / 2
    else
        r = sqrt(1 + 4*x*y)
        return x + y - r + log((1 + r)/(2*y))
    end
end

# ζ from eq. (56)
ζ(x, y) = copysign(sqrt(2*half_ζ2(x,y)), x + 1 - y)

# θ/sin(θ) and derivative·sin(θ) for small θ
theta_over_sin(θ) = θ < 1e-4 ? 1 + θ^2/6 : θ/sin(θ)
theta_prime_sin(θ) = θ < 1e-4 ? θ^2*(θ^2/45 + 1/3) : 1 - θ/tan(θ)

# ρ from eq. (98)
ρ(θo, ξ) = sqrt(θo^2 + ξ^2)

# r and r′·sin(θ) from eq. (96)
function r(θ, y, ξ)
    tos = theta_over_sin(θ); xs = ξ/tos
    (1 + sqrt(1 + xs^2))*tos/(2*y)
end

function r_prime_sin(θ, y, ξ)
    tos = theta_over_sin(θ); xs = ξ/tos
    (1 + 1/sqrt(1 + xs^2))*theta_prime_sin(θ)/(2*y)
end

# f from eq. (96)
f(θ, y, ξ) = begin
    r0 = r(θ, y, ξ); d = r0 - cos(θ)
    (r_prime_sin(θ, y, ξ) - d*r0)/(d^2 + sin(θ)^2)
end

# ψ from eq. (97)
ψ(θ, ξ) = begin
    tos = theta_over_sin(θ); rv = ρ(tos, ξ)
    cos(θ)*rv - sqrt(1 + ξ^2) - log((tos + rv)/(1 + sqrt(1 + ξ^2)))
end

# f1, f2 from eq. (100)
f1_f2(x, M) = (x + M - sqrt(4*x + 2*M), x + M + sqrt(4*x + 2*M))

# series expansion (section 3)
function MarcumQ_small_x(M, x, y)
    s = 0.0; tf = 1.0; n = 0
    while true
        t = tf*Q(M+n, y); s += t
        abs(t) <= 1e-17*abs(s) && break
        n += 1; tf *= x/n
    end
    exp(-x)*s
end

# asymptotic expansion for large x·y (section 4.1)
function MarcumQ_large_xy(M, x, y, ξ)
    δ = sqrt(y) - sqrt(x); σ = δ^2/ξ; ρ0 = sqrt(y/x)
    ρfac = (y/x)^(M/2)/sqrt(8*π); ef = exp(-δ^2)*sqrt(ξ)
    Φ = abs(δ)<1e-5 ? (σ==0 ? 0.0 : sqrt(π/σ)-2*sqrt(ξ)) : sqrt(π/σ)*erfc(abs(δ))
    Ψ = ρ0==1 ? 0.5 : copysign(ρ0^(M-0.5)/2*erfc(abs(δ)), ρ0-1)
    s = x>y ? 1 : 0; n = 0; ρt = ρfac
    while true
        s += Ψ; abs(Ψ) <= 1e-17*abs(s) && break
        n += 1; ρt = -ρt; ef /= ξ
        Φ = (ef - σ*Φ)/(n - 0.5)
        lnAn = lnA(n, M-1)
        Ψ = ρt*exp(lnAn)*(1 - exp(lnA(n,M)-lnAn)/ρ0)*Φ
    end
    max(s, 0)
end

# recurrence (eq. 14)
function MarcumQ_recurrence(M, x, y, ξ)
    root = sqrt(y/x)
    μ = M - ceil(M - (sqrt(2*ξ) - 1))
    Qm1 = MarcumQ_modified(μ-1, x, y); Q0 = MarcumQ_modified(μ, x, y)
    while μ < M - 1e-15*M
        cm = root*c_μ(μ, ξ)
        Q1 = (1 + cm)*Q0 - cm*Qm1
        Qm1, Q0 = Q0, Q1; μ += 1
    end
    Q0
end

# asymptotic for large M (section 4.2)
function MarcumQ_large_M(M, x, y)
    ζv = ζ(x, y); ehalf = exp(-M*half_ζ2(x,y))
    Ψ = zeros(Float64, 20)
    Ψ[1] = sqrt(π/(2*M))*erfc(-ζv*sqrt(M/2)); Ψ[2] = ehalf/M
    s = 0.0; k = 1
    while true
        Bk = sum(Ψ[j]/M^(k-j) for j in 1:k); s += Bk
        abs(Bk) <= 1e-17*abs(s) && break
        k += 1; Ψ[k] = (k-1)/M*Ψ[k-2] + (-ζv)^(k-1)/M*ehalf
    end
    erfc(-ζv*sqrt(M/2))/2 - sqrt(M/(2*π))*s
end

# quadrature (section 5)
function MarcumQ_quadrature(M, x, y, ξ)
    I, _ = quadgk(θ->exp(M*ψ(θ,ξ))*f(θ,y,ξ), 0, π-1/512; rtol=1e-12)
    I *= exp(-M*half_ζ2(x,y))/π
    x + 1 < y ? I : 1 + I
end

# core modified Marcum Q
function MarcumQ_modified(M, x, y)
    @assert 1 ≤ M ≤ 1e4; @assert 0 ≤ x ≤ 1e4; @assert 0 ≤ y ≤ 1e4
    ξ = 2*sqrt(x*y); f1, f2 = f1_f2(x, M)
    Qv = x < 30 ? MarcumQ_small_x(M,x,y) :
        (ξ>30 && M^2<2*ξ) ? MarcumQ_large_xy(M,x,y,ξ) :
        (f1<y<f2 && M<135)  ? MarcumQ_recurrence(M,x,y,ξ) :
        (f1<y<f2 && M≥135)  ? error("Large-M not implemented") :
                               MarcumQ_quadrature(M,x/M,y/M,ξ/M)
    Qv > 1 && Qv < 1+1e-12 && (Qv = 1)
    @assert 0 ≤ Qv ≤ 1; Qv
end

"""
    MarcumQ(μ::Float64, a::Float64, b::Float64)

Compute the generalized Marcum Q-function of order `μ` with non-centrality parameter `a` and threshold `b`.

Reference:
    [1] https://arxiv.org/pdf/1311.0681v1
"""
function MarcumQ(μ::Float64, a::Float64, b::Float64)
    return MarcumQ_modified(μ, a^2/2.0, b^2/2.0)
end

function MarcumQ(μ::Real, a::Real, b::Real)
    return MarcumQ_modified(Float64(μ), Float64(a)^2/2.0, Float64(b)^2/2.0)
end

function MarcumQ(μ::Real, a::Real, b::AbstractArray{<:Real})
    return [MarcumQ(μ, a, bᵢ) for bᵢ in b]
end

function MarcumQ(μ::Real, a::AbstractArray{<:Real}, b::Real)
    return [MarcumQ(μ, aᵢ, b) for aᵢ in a]
end

function MarcumQ(μ::AbstractArray{<:Real}, a::Real, b::Real)
    return [MarcumQ(μᵢ, a, b) for μᵢ in μ]
end

function MarcumQ(a::Real, b::Real)
    return MarcumQ(1, a, b)
end
