using SpecialFunctions: gamma

export U, V, W, dU, dV, dW

"""
    U(a::T, x::T) where {T <: AbstractFloat}

Compute the parabolic cylinder function U(a,x) of the first kind for real parameters.
Supports any `AbstractFloat` type (e.g., `Float32`, `Float64`, `BigFloat`).

S. Zhang and J. Jin, 'Computation of Special functions' (Wiley, 1966), E. Cojocaru, January 2009
"""
function U(a::T, x::T) where {T <: AbstractFloat}
    ε = eps(T) * 10

    return if abs(x) ≤ 5
        c = zeros(T, 100)
        c[1] = a
        c₀, c₁ = one(T), a
        for k in 4:2:200
            m = k ÷ 2
            c[m] = a * c₁ + T(k - 2) * T(k - 3) * c₀ / 4
            c₀, c₁ = c₁, c[m]
        end

        y₁ = one(T)
        term = one(T)
        @inbounds for k in 1:100
            term *= x^2 / (T(2) * T(k) * T(2k - 1))
            Δ = c[k] * term
            y₁ += Δ
            if abs(Δ / y₁) ≤ ε && k > 30
                break
            end
        end

        d = zeros(T, 100)
        d[1], d[2] = one(T), a
        d₁, d₂ = one(T), a
        for k in 5:2:160
            m = (k + 1) ÷ 2
            d[m] = a * d₂ + T(k - 2) * T(k - 3) * d₁ / 4
            d₁, d₂ = d₂, d[m]
        end

        y₂ = one(T)
        term = one(T)
        @inbounds for k in 1:100
            term *= x^2 / (T(2) * T(k) * T(2k + 1))
            Δ = d[k + 1] * term
            y₂ += Δ
            if abs(Δ / y₂) ≤ ε && k > 30
                break
            end
        end
        y₂ *= x

        if a < 0 && isapprox(a + T(0.5), round(a + T(0.5)))
            θ = T(π) * (T(0.25) + a / 2)
            f₁ = gamma(T(0.25) - a / 2) / (sqrt(T(π)) * T(2)^(a / 2 + T(0.25)))
            f₂ = gamma(T(0.75) - a / 2) / (sqrt(T(π)) * T(2)^(a / 2 - T(0.25)))
            return cos(θ) * f₁ * y₁ - sin(θ) * f₂ * y₂
        else
            prefactor = sqrt(T(π)) / T(2)^(a / 2 + T(0.25))
            g₁ = gamma(T(0.25) + a / 2)
            g₃ = gamma(T(0.75) + a / 2)
            return prefactor * (y₁ / g₃ - sqrt(T(2)) * y₂ / g₁)
        end
    else
        q = exp(-x^2 / 4)
        a₀ = q * abs(x)^(-a - T(0.5))
        r = one(T)
        u = one(T)
        for k in 1:20
            r *= -(T(2k) + a - T(0.5)) * (T(2k) + a - T(1.5)) / (T(2k) * x^2)
            u += r
            if abs(r / u) < ε
                break
            end
        end
        u *= a₀

        if x < 0
            if a < 0 && isapprox(a + T(0.5), round(a + T(0.5)))
                return -u * sinpi(a)
            else
                v = V(a, -x)
                return T(π) * v / gamma(a + T(0.5)) - sinpi(a) * u
            end
        else
            return u
        end
    end
end

"""
    V(a::T, x::T) where {T <: AbstractFloat}

Compute the parabolic cylinder function V(a,x).
Supports any `AbstractFloat` type (e.g., `Float32`, `Float64`, `BigFloat`).
"""
function V(a::T, x::T) where {T <: AbstractFloat}
    ε = eps(T) * 10

    c = zeros(T, 100)
    c[1] = a
    c₀, c₁ = one(T), a
    for k in 4:2:200
        m = k ÷ 2
        c[m] = a * c₁ + T(k - 2) * T(k - 3) * c₀ / 4
        c₀, c₁ = c₁, c[m]
    end

    y₁ = one(T)
    term = one(T)
    @inbounds for k in 1:100
        term *= x^2 / (T(2) * T(k) * T(2k - 1))
        Δ = c[k] * term
        y₁ += Δ
        if abs(Δ / y₁) ≤ ε && k > 30
            break
        end
    end

    d = zeros(T, 100)
    d[1], d[2] = one(T), a
    d₁, d₂ = one(T), a
    for k in 5:2:160
        m = (k + 1) ÷ 2
        d[m] = a * d₂ + T(k - 2) * T(k - 3) * d₁ / 4
        d₁, d₂ = d₂, d[m]
    end

    y₂ = one(T)
    term = one(T)
    @inbounds for k in 1:100
        term *= x^2 / (T(2) * T(k) * T(2k + 1))
        Δ = d[k + 1] * term
        y₂ += Δ
        if abs(Δ / y₂) ≤ ε && k > 30
            break
        end
    end
    y₂ *= x

    if a < 0 && isapprox(a + T(0.5), round(a + T(0.5)))
        θ = T(π) * (T(0.25) + a / 2)
        f₁ = gamma(T(0.25) - a / 2) / (sqrt(T(π)) * T(2)^(a / 2 + T(0.25)))
        f₂ = gamma(T(0.75) - a / 2) / (sqrt(T(π)) * T(2)^(a / 2 - T(0.25)))
        v = (sin(θ) * f₁ * y₁ + cos(θ) * f₂ * y₂) / gamma(T(0.5) - a)
    else
        sinπa = sinpi(a)
        g₀ = gamma(T(0.5) + a)
        p₀ = g₀ / (sqrt(T(π)) * T(2)^(a / 2 + T(0.25)))
        g₁ = gamma(T(0.25) + a / 2)
        g₃ = gamma(T(0.75) + a / 2)
        v = p₀ * (y₁ * (one(T) + sinπa) / g₃ + sqrt(T(2)) * y₂ * (one(T) - sinπa) / g₁)
    end

    return v
end

"""
    W(a::T, x::T) where {T <: AbstractFloat}

Compute the parabolic cylinder function W(a,x) for real parameters.
Supports any `AbstractFloat` type (e.g., `Float32`, `Float64`, `BigFloat`).
"""
function W(a::T, x::T) where {T <: AbstractFloat}
    ε = eps(T) * 10
    if abs(x) ≤ 8
        p₀ = T(2)^(-T(3) / 4)
        g₁ = abs(gamma(Complex{T}(T(1) / 4, a / 2)))
        g₃ = abs(gamma(Complex{T}(T(3) / 4, a / 2)))
        f₁ = sqrt(g₁ / g₃)
        f₂ = sqrt(T(2) * g₃ / g₁)

        c = zeros(T, 100)
        c[1] = a
        c₀, c₁ = one(T), a
        for k in 4:2:200
            m = k ÷ 2
            c[m] = a * c₁ - T(k - 2) * T(k - 3) * c₀ / 4
            c₀, c₁ = c₁, c[m]
        end

        y₁ = one(T)
        term = one(T)
        for k in 1:100
            term *= x^2 / (T(2) * T(k) * T(2k - 1))
            Δ = c[k] * term
            y₁ += Δ
            if abs(Δ / y₁) ≤ ε && k > 30
                break
            end
        end

        d = zeros(T, 100)
        d[1], d[2] = one(T), a
        d₁, d₂ = one(T), a
        for k in 5:2:160
            m = (k + 1) ÷ 2
            d[m] = a * d₂ - T(k - 2) * T(k - 3) * d₁ / 4
            d₁, d₂ = d₂, d[m]
        end

        y₂ = one(T)
        term = one(T)
        for k in 1:100
            term *= x^2 / (T(2) * T(k) * T(2k + 1))
            Δ = d[k + 1] * term
            y₂ += Δ
            if abs(Δ / y₂) ≤ ε && k > 30
                break
            end
        end
        y₂ *= x

        return p₀ * (f₁ * y₁ - f₂ * y₂)

    else
        u = zeros(T, 21)
        v = zeros(T, 21)

        g₀ = gamma(Complex{T}(T(1) / 2, a))
        ϕ₂ = imag(g₀)

        gref = gamma(Complex{T}(T(1) / 2, a))
        gr₀, gi₀ = real(gref), imag(gref)
        den = gr₀^2 + gi₀^2

        for k in 2:2:40
            m = k ÷ 2
            g = gamma(Complex{T}(T(k) + T(0.5), a))
            gr, gi = real(g), imag(g)
            u[m] = (gr * gr₀ + gi * gi₀) / den
            v[m] = (gr₀ * gi - gr * gi₀) / den
        end

        x² = x^2
        sv₁ = v[1] / (T(2) * x²)
        su₂ = u[1] / (T(2) * x²)
        fac = one(T)
        for k in 3:2:19
            fac *= -T(k) * T(k - 1)
            denom = fac * T(2)^k * x^(2k)
            sv₁ += v[k] / denom
            su₂ += u[k] / denom
        end

        sv₂ = zero(T)
        su₁ = zero(T)
        fac = one(T)
        for k in 2:2:20
            fac *= -T(k) * T(k - 1)
            denom = fac * T(2)^k * x^(2k)
            sv₂ += v[k] / denom
            su₁ += u[k] / denom
        end

        s₁ = one(T) + sv₁ - su₁
        s₂ = -sv₂ - su₂

        ea = exp(T(π) * a)
        S = sqrt(one(T) + ea^2)
        f₊ = S + ea
        f₋ = S - ea

        ϕ = x^2 / 4 - a * log(abs(x)) + T(π) / 4 + ϕ₂ / 2

        if x > 0
            return sqrt(T(2) * f₋ / abs(x)) * (s₁ * cos(ϕ) - s₂ * sin(ϕ))
        else
            return sqrt(T(2) * f₊ / abs(x)) * (s₁ * sin(ϕ) + s₂ * cos(ϕ))
        end
    end
end

"""
    dU(a::T, x::T) where {T <: AbstractFloat}

Compute the derivative of the parabolic cylinder function U(a,x) for real parameters.
Supports any `AbstractFloat` type (e.g., `Float32`, `Float64`, `BigFloat`).
"""
function dU(a::T, x::T) where {T <: AbstractFloat}
    ε = eps(T) * 10

    c = zeros(T, 101)
    c[1] = a
    c₀, c₁ = one(T), a
    for k in 4:2:202
        m = k ÷ 2
        c[m] = a * c₁ + T(k - 2) * T(k - 3) * c₀ / 4
        c₀, c₁ = c₁, c[m]
    end

    y₁ = a
    term = one(T)
    @inbounds for k in 1:100
        term *= x^2 / (T(2) * T(k) * T(2k + 1))
        Δ = c[k + 1] * term
        y₁ += Δ
        if abs(Δ / y₁) ≤ ε && k > 30
            break
        end
    end
    y₁ *= x

    d = zeros(T, 101)
    d[1], d[2] = one(T), a
    d₁, d₂ = one(T), a
    for k in 5:2:160
        m = (k + 1) ÷ 2
        d[m] = a * d₂ + T(k - 2) * T(k - 3) * d₁ / 4
        d₁, d₂ = d₂, d[m]
    end

    y₂ = one(T)
    term = one(T)
    @inbounds for k in 1:100
        term *= x^2 / (T(2) * T(k) * T(2k - 1))
        Δ = d[k + 1] * term
        y₂ += Δ
        if abs(Δ / y₂) ≤ ε && k > 30
            break
        end
    end

    if a < 0 && isapprox(a + T(0.5), round(a + T(0.5)))
        θ = T(π) * (T(0.25) + a / 2)
        f₁ = gamma(T(0.25) - a / 2) / (sqrt(T(π)) * T(2)^(a / 2 + T(0.25)))
        f₂ = gamma(T(0.75) - a / 2) / (sqrt(T(π)) * T(2)^(a / 2 - T(0.25)))
        du = cos(θ) * f₁ * y₁ - sin(θ) * f₂ * y₂
    else
        prefactor = sqrt(T(π)) / T(2)^(a / 2 + T(0.25))
        g₁ = gamma(T(0.25) + a / 2)
        g₃ = gamma(T(0.75) + a / 2)
        du = prefactor * (y₁ / g₃ - sqrt(T(2)) * y₂ / g₁)
    end

    return du
end

"""
    dV(a::T, x::T) where {T <: AbstractFloat}

Compute the derivative of the parabolic cylinder function V(a,x) for real parameters.
Supports any `AbstractFloat` type (e.g., `Float32`, `Float64`, `BigFloat`).
"""
function dV(a::T, x::T) where {T <: AbstractFloat}
    ε = eps(T) * 10

    c = zeros(T, 101)
    c[1] = a
    c₀, c₁ = one(T), a
    for k in 4:2:202
        m = k ÷ 2
        c[m] = a * c₁ + T(k - 2) * T(k - 3) * c₀ / 4
        c₀, c₁ = c₁, c[m]
    end

    y₁ = a
    term = one(T)
    @inbounds for k in 1:100
        term *= x^2 / (T(2) * T(k) * T(2k + 1))
        Δ = c[k + 1] * term
        y₁ += Δ
        if abs(Δ / y₁) ≤ ε && k > 30
            break
        end
    end
    y₁ *= x

    d = zeros(T, 101)
    d[1], d[2] = one(T), a
    d₁, d₂ = one(T), a
    for k in 5:2:160
        m = (k + 1) ÷ 2
        d[m] = a * d₂ + T(k - 2) * T(k - 3) * d₁ / 4
        d₁, d₂ = d₂, d[m]
    end

    y₂ = one(T)
    term = one(T)
    @inbounds for k in 1:100
        term *= x^2 / (T(2) * T(k) * T(2k - 1))
        Δ = d[k + 1] * term
        y₂ += Δ
        if abs(Δ / y₂) ≤ ε && k > 30
            break
        end
    end

    if a < 0 && isapprox(a + T(0.5), round(a + T(0.5)))
        θ = T(π) * (T(0.25) + a / 2)
        f₁ = gamma(T(0.25) - a / 2) / (sqrt(T(π)) * T(2)^(a / 2 + T(0.25)))
        f₂ = gamma(T(0.75) - a / 2) / (sqrt(T(π)) * T(2)^(a / 2 - T(0.25)))
        dv = (sin(θ) * f₁ * y₁ + cos(θ) * f₂ * y₂) / gamma(T(0.5) - a)
    else
        sinπa = sinpi(a)
        g₀ = gamma(T(0.5) + a)
        p₀ = g₀ / (sqrt(T(π)) * T(2)^(a / 2 + T(0.25)))
        g₁ = gamma(T(0.25) + a / 2)
        g₃ = gamma(T(0.75) + a / 2)
        dv = p₀ * (y₁ * (one(T) + sinπa) / g₃ + sqrt(T(2)) * y₂ * (one(T) - sinπa) / g₁)
    end

    return dv
end


"""
    dW(a::T, x::T) where {T <: AbstractFloat}

Compute the derivative of the parabolic cylinder function W with parameters `a` evaluated at `x`.
Supports any `AbstractFloat` type (e.g., `Float32`, `Float64`, `BigFloat`).
"""
function dW(a::T, x::T) where {T <: AbstractFloat}
    ε = eps(T) * 10
    p₀ = T(2)^(-T(3) / 4)

    g₁ = abs(gamma(Complex{T}(T(1) / 4, a / 2)))
    g₃ = abs(gamma(Complex{T}(T(3) / 4, a / 2)))
    f₁ = sqrt(g₁ / g₃)
    f₂ = sqrt(T(2) * g₃ / g₁)

    c = zeros(T, 101)
    c[1] = a
    c₀, c₁ = one(T), a
    for k in 4:2:202
        m = k ÷ 2
        c[m] = a * c₁ - T(k - 2) * T(k - 3) * c₀ / 4
        c₀, c₁ = c₁, c[m]
    end

    y₁ = a
    term = one(T)
    @inbounds for k in 1:100
        term *= x^2 / (T(2) * T(k) * T(2k + 1))
        Δ = c[k + 1] * term
        y₁ += Δ
        if abs(Δ / y₁) ≤ ε && k > 30
            break
        end
    end
    y₁ *= x

    d = zeros(T, 101)
    d[1], d[2] = one(T), a
    d₁, d₂ = one(T), a
    for k in 5:2:160
        m = (k + 1) ÷ 2
        d[m] = a * d₂ - T(k - 2) * T(k - 3) * d₁ / 4
        d₁, d₂ = d₂, d[m]
    end

    y₂ = one(T)
    term = one(T)
    @inbounds for k in 1:100
        term *= x^2 / (T(2) * T(k) * T(2k - 1))
        Δ = d[k + 1] * term
        y₂ += Δ
        if abs(Δ / y₂) ≤ ε && k > 30
            break
        end
    end

    return p₀ * (f₁ * y₁ - f₂ * y₂)
end

# Promotion methods for Real inputs
for func in (:U, :V, :W, :dU, :dV, :dW)
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
