using SpecialFunctions: besseli, gamma_inc, loggamma, erfc

export MarcumQ, dQdb

Q(μ::T, x::T) where {T<:Number} = gamma_inc(μ, x)[2]

# continued‐fraction ratio for scaled I_μ
c_μ(μ::T, ξ::T) where {T<:Number} = besseli(μ, ξ) / besseli(μ-1, ξ)

# log of Aₙ from eq. (32)
lnA(n::Integer, μ::T) where {T<:Number} = 
    loggamma(μ + T(0.5) + n) - loggamma(μ + T(0.5) - n) - n*log(T(2)) - loggamma(n+1)

# ζ²/2 from eq. (84)
function half_ζ2(x::T, y::T) where {T<:Number}
    δ = y - x - one(T)
    if abs(δ) < T(1e-3)
        c = T[1,
             -(T(3)*x + one(T))/T(3),
             ((T(72)*x + T(42))*x + T(7))/T(36),
             -(((T(2700)*x + T(2142))*x + T(657))*x + T(73))/T(540),
             ((((T(181440)*x + T(177552))*x + T(76356))*x + T(15972))*x + T(1331))/T(12960)]
        z = δ / ( (T(2)*x + one(T))^2 )
        s = sum(c[k] * z^k for k in eachindex(c))
        return (T(2)*x + one(T))^3 * s^2 / T(2)
    else
        r = sqrt(one(T) + T(4)*x*y)
        return x + y - r + log((one(T) + r)/(T(2)*y))
    end
end

# ζ from eq. (56)
ζ(x::T, y::T) where {T<:Number} = copysign(sqrt(T(2)*half_ζ2(x,y)), x + one(T) - y)

# θ/sin(θ) and derivative·sin(θ) for small θ
theta_over_sin(θ::T) where {T<:Number} = θ < T(1e-4) ? one(T) + θ^2/T(6) : θ/sin(θ)
theta_prime_sin(θ::T) where {T<:Number} = θ < T(1e-4) ? θ^2*(θ^2/T(45) + one(T)/T(3)) : one(T) - θ/tan(θ)

# ρ from eq. (98)
ρ(θo::T, ξ::T) where {T<:Number} = sqrt(θo^2 + ξ^2)

# r and r′·sin(θ) from eq. (96)
function r(θ::T, y::T, ξ::T) where {T<:Number}
    tos = theta_over_sin(θ); xs = ξ/tos
    (one(T) + sqrt(one(T) + xs^2))*tos/(T(2)*y)
end

function r_prime_sin(θ::T, y::T, ξ::T) where {T<:Number}
    tos = theta_over_sin(θ); xs = ξ/tos
    (one(T) + one(T)/sqrt(one(T) + xs^2))*theta_prime_sin(θ)/(T(2)*y)
end

# f from eq. (96)
f(θ::T, y::T, ξ::T) where {T<:Number} = begin
    r0 = r(θ, y, ξ); d = r0 - cos(θ)
    (r_prime_sin(θ, y, ξ) - d*r0)/(d^2 + sin(θ)^2)
end

# ψ from eq. (97)
ψ(θ::T, ξ::T) where {T<:Number} = begin
    tos = theta_over_sin(θ); rv = ρ(tos, ξ)
    cos(θ)*rv - sqrt(one(T) + ξ^2) - log((tos + rv)/(one(T) + sqrt(one(T) + ξ^2)))
end

# f1, f2 from eq. (100)
f1_f2(x::T, M::T) where {T<:Number} = (x + M - sqrt(T(4)*x + T(2)*M), x + M + sqrt(T(4)*x + T(2)*M))

# series expansion (section 3)
function MarcumQ_small_x(M::T, x::T, y::T) where {T<:Number}
    s = zero(T); tf = one(T); n = 0
    while true
        t = tf * Q(M + n, y); s += t
        abs(t) <= eps(T)*abs(s) && break
        n += 1; tf *= x/T(n)
    end
    exp(-x)*s
end

# asymptotic expansion for large x·y (section 4.1)
function MarcumQ_large_xy(M::T, x::T, y::T, ξ::T) where {T<:Number}
    δ = sqrt(y) - sqrt(x); σ = δ^2/ξ; ρ0 = sqrt(y/x)
    ρfac = (y/x)^(M/T(2))/sqrt(T(8)*π); ef = exp(-δ^2)*sqrt(ξ)
    Φ = abs(δ)<T(1e-5) ? (σ==zero(T) ? zero(T) : sqrt(π/σ)-T(2)*sqrt(ξ)) : sqrt(π/σ)*erfc(abs(δ))
    Ψ = ρ0==one(T) ? half(T) : copysign(ρ0^(M-T(0.5))/T(2)*erfc(abs(δ)), ρ0-one(T))
    s = x>y ? one(T) : zero(T); n = 0; ρt = ρfac
    while true
        s += Ψ; abs(Ψ) <= eps(T)*abs(s) && break
        n += 1; ρt = -ρt; ef /= ξ
        Φ = (ef - σ*Φ)/(n - T(0.5))
        lnAn = lnA(n, M - one(T))
        Ψ = ρt*exp(lnAn)*(one(T) - exp(lnA(n,M)-lnAn)/ρ0)*Φ
    end
    max(s, zero(T))
end

# recurrence (eq. 14)
function MarcumQ_recurrence(M::T, x::T, y::T, ξ::T) where {T<:Number}
    root = sqrt(y/x)
    μ = M - ceil(M - (sqrt(T(2)*ξ) - one(T)))
    Qm1 = MarcumQ_modified(μ-one(T), x, y); Q0 = MarcumQ_modified(μ, x, y)
    while μ < M - eps(T)*M
        cm = root*c_μ(μ, ξ)
        Q1 = (one(T)+cm)*Q0 - cm*Qm1
        Qm1, Q0 = Q0, Q1; μ += one(T)
    end
    Q0
end

# asymptotic for large M (section 4.2)
function MarcumQ_large_M(M::T, x::T, y::T) where {T<:Number}
    ζv = ζ(x, y); ehalf = exp(-M*half_ζ2(x,y))
    Ψ = [zero(T) for _ in 1:20]
    Ψ[1] = sqrt(π/(T(2)*M))*erfc(-ζv*sqrt(M/T(2))); Ψ[2] = ehalf/M
    s = zero(T); k = one(Int)
    while true
        Bk = sum(Ψ[j]/M^(k-j) for j in one(Int):k); s += Bk
        abs(Bk) <= eps(T)*abs(s) && break
        k += one(Int); Ψ[k] = (k-one(Int))/M*Ψ[k-one(Int)] + (-ζv)^(k-one(Int))/M*ehalf
    end
    erfc(-ζv*sqrt(M/T(2)))/T(2) - sqrt(M/(T(2)*π))*s
end

# quadrature (section 5)
function MarcumQ_quadrature(M::T, x::T, y::T, ξ::T) where {T<:Number}
    I, _ = quadgk(θ->exp(M*ψ(θ,ξ))*f(θ,y,ξ), zero(T), π-one(T)/T(512); rtol=eps(T))
    I *= exp(-M*half_ζ2(x,y))/π
    x + one(T) < y ? I : one(T) + I
end

# core modified Marcum Q
function MarcumQ_modified(M::T, x::T, y::T) where {T<:Number}
    ξ = T(2)*sqrt(x*y); f1, f2 = f1_f2(x, M)
    Qv = x < T(30) ? MarcumQ_small_x(M,x,y) :
         (ξ>T(30) && M^2<T(2)*ξ) ? MarcumQ_large_xy(M,x,y,ξ) :
         (f1<y<f2 && M<T(135))   ? MarcumQ_recurrence(M,x,y,ξ) :
         (f1<y<f2 && M>=T(135))  ? error("Large-M not implemented") :
                                   MarcumQ_quadrature(M,x/T(M),y/T(M),ξ/T(M))
    Qv > one(T) && Qv < one(T)+eps(T) && (Qv = one(T))
    Qv
end

# standard Marcum Q
MarcumQ(M::T, a::T, b::T) where {T<:Number} = MarcumQ_modified(M, a^2/T(2), b^2/T(2))

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

"""
    dQdb(M, a, b)

Derivative ∂Q_M(a,b)/∂b of the (standard) Marcum Q-function of order `M`.
Requires `M` integer ≥1 and `a>0`.
"""
function dQdb(M::Integer, a::T, b::T) where {T<:Number}
    @assert M ≥ 1 "order M must be ≥1"
    @assert a != zero(a) "a must be nonzero"
    coeff = b^M / a^(M-1)
    return -coeff * exp( -(a^2 + b^2)/T(2) ) * besseli(M-1, a*b)
end

function dQdb(M::Real, a::Real, b::Real)
    return dQdb(Float64(M), Float64(a), Float64(b))
end

# Array handling for b parameter
function dQdb(M::Real, a::Real, b::AbstractArray{<:Real})
    return [dQdb(M, a, bᵢ) for bᵢ in b]
end

# Array handling for a parameter
function dQdb(M::Real, a::AbstractArray{<:Real}, b::Real)
    return [dQdb(M, aᵢ, b) for aᵢ in a]
end

# Array handling for M parameter
function dQdb(M::AbstractArray{<:Real}, a::Real, b::Real)
    return [dQdb(Mᵢ, a, b) for Mᵢ in M]
end

# Convenience method for M=1
function dQdb(a::Real, b::Real)
    return dQdb(1, a, b)
end