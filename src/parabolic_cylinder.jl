using SpecialFunctions: gamma

export U, V, W, dU, dV, dW

"""
    U(a::Float64, x::Float64)::Float64

Compute the parabolic cylinder function U(a,x) of the first kind for real parameters.
    
S. Zhang and J. Jin, 'Computation of Special functions' (Wiley, 1966), E. Cojocaru, January 2009
"""
function U(a::Float64, x::Float64)::Float64
    ε = 1e-15

    if abs(x) ≤ 5
        c = zeros(Float64, 100)
        c[1] = a
        c₀, c₁ = 1.0, a
        for k in 4:2:200
            m = k ÷ 2
            c[m] = a * c₁ + (k - 2) * (k - 3) * c₀ / 4
            c₀, c₁ = c₁, c[m]
        end

        y₁ = 1.0
        term = 1.0
        @inbounds for k in 1:100
            term *= 0.5 * x^2 / (k * (2k - 1))
            Δ = c[k] * term
            y₁ += Δ
            if abs(Δ / y₁) ≤ ε && k > 30
                break
            end
        end

        d = zeros(Float64, 100)
        d[1], d[2] = 1.0, a
        d₁, d₂ = 1.0, a
        for k in 5:2:160
            m = (k + 1) ÷ 2
            d[m] = a * d₂ + 0.25 * (k - 2) * (k - 3) * d₁
            d₁, d₂ = d₂, d[m]
        end

        y₂ = 1.0
        term = 1.0
        @inbounds for k in 1:100
            term *= 0.5 * x^2 / (k * (2k + 1))
            Δ = d[k + 1] * term
            y₂ += Δ
            if abs(Δ / y₂) ≤ ε && k > 30
                break
            end
        end
        y₂ *= x

        if a < 0 && isapprox(a + 0.5, round(a + 0.5))
            θ = pi * (0.25 + a / 2)
            f₁ = gamma(0.25 - a / 2) / (√π * 2^(a / 2 + 0.25))
            f₂ = gamma(0.75 - a / 2) / (√π * 2^(a / 2 - 0.25))
            return cos(θ) * f₁ * y₁ - sin(θ) * f₂ * y₂
        else
            prefactor = √π / 2^(a / 2 + 0.25)
            g₁ = gamma(0.25 + a / 2)
            g₃ = gamma(0.75 + a / 2)
            return prefactor * (y₁ / g₃ - √2 * y₂ / g₁)
        end
    else
        q = exp(-x^2 / 4)
        a₀ = q * abs(x)^(-a - 0.5)
        r = 1.0
        u = 1.0
        for k in 1:20
            r *= -(2k + a - 0.5) * (2k + a - 1.5) / (2k * x^2)
            u += r
            if abs(r / u) < ε
                break
            end
        end
        u *= a₀

        if x < 0
            if a < 0 && isapprox(a + 0.5, round(a + 0.5))
                return -u * sinpi(a)
            else
                v = V(a, -x)
                return pi * v / gamma(a + 0.5) - sinpi(a) * u
            end
        else
            return u
        end
    end
end

"""
    V(a::Float64, x::Float64)::Float64

Compute the parabolic cylinder function V(a,x).
"""
function V(a::Float64, x::Float64)::Float64
    ε = 1e-15

    c = zeros(Float64, 100)
    c[1] = a
    c₀, c₁ = 1.0, a
    for k in 4:2:200
        m = k ÷ 2
        c[m] = a * c₁ + (k - 2) * (k - 3) * c₀ / 4
        c₀, c₁ = c₁, c[m]
    end

    y₁ = 1.0
    term = 1.0
    @inbounds for k in 1:100
        term *= 0.5 * x^2 / (k * (2k - 1))
        Δ = c[k] * term
        y₁ += Δ
        if abs(Δ / y₁) ≤ ε && k > 30
            break
        end
    end

    d = zeros(Float64, 100)
    d[1], d[2] = 1.0, a
    d₁, d₂ = 1.0, a
    for k in 5:2:160
        m = (k + 1) ÷ 2
        d[m] = a * d₂ + 0.25 * (k - 2) * (k - 3) * d₁
        d₁, d₂ = d₂, d[m]
    end

    y₂ = 1.0
    term = 1.0
    @inbounds for k in 1:100
        term *= 0.5 * x^2 / (k * (2k + 1))
        Δ = d[k + 1] * term
        y₂ += Δ
        if abs(Δ / y₂) ≤ ε && k > 30
            break
        end
    end
    y₂ *= x

    if a < 0 && isapprox(a + 0.5, round(a + 0.5))
        θ = pi * (0.25 + a / 2)
        f₁ = gamma(0.25 - a / 2) / (√π * 2^(a / 2 + 0.25))
        f₂ = gamma(0.75 - a / 2) / (√π * 2^(a / 2 - 0.25))
        v = (sin(θ) * f₁ * y₁ + cos(θ) * f₂ * y₂) / gamma(0.5 - a)
    else
        sinπa = sinpi(a)
        g₀ = gamma(0.5 + a)
        p₀ = g₀ / (√π * 2^(a / 2 + 0.25))
        g₁ = gamma(0.25 + a / 2)
        g₃ = gamma(0.75 + a / 2)
        v = p₀ * (y₁ * (1 + sinπa) / g₃ + √2 * y₂ * (1 - sinπa) / g₁)
    end

    return v
end

"""
    W(a::Float64, x::Float64)::Float64

Compute the parabolic cylinder function W(a,x) for real parameters.
"""
function W(a::Float64, x::Float64)::Float64
    ε = 1e-15
    if abs(x) ≤ 8
        p₀ = 2^(-3/4)
        g₁ = abs(gamma(Complex(1/4, a/2)))
        g₃ = abs(gamma(Complex(3/4, a/2)))
        f₁ = sqrt(g₁ / g₃)
        f₂ = sqrt(2 * g₃ / g₁)

        c = zeros(Float64, 100)
        c[1] = a
        c₀, c₁ = 1.0, a
        for k in 4:2:200
            m = k ÷ 2
            c[m] = a * c₁ - (k - 2) * (k - 3) * c₀ / 4
            c₀, c₁ = c₁, c[m]
        end

        y₁ = 1.0
        term = 1.0
        for k in 1:100
            term *= 0.5 * x^2 / (k * (2k - 1))
            Δ = c[k] * term
            y₁ += Δ
            if abs(Δ / y₁) ≤ ε && k > 30
                break
            end
        end

        d = zeros(Float64, 100)
        d[1], d[2] = 1.0, a
        d₁, d₂ = 1.0, a
        for k in 5:2:160
            m = (k + 1) ÷ 2
            d[m] = a * d₂ - 0.25 * (k - 2) * (k - 3) * d₁
            d₁, d₂ = d₂, d[m]
        end

        y₂ = 1.0
        term = 1.0
        for k in 1:100
            term *= 0.5 * x^2 / (k * (2k + 1))
            Δ = d[k + 1] * term
            y₂ += Δ
            if abs(Δ / y₂) ≤ ε && k > 30
                break
            end
        end
        y₂ *= x

        return p₀ * (f₁ * y₁ - f₂ * y₂)

    else
        u = zeros(Float64, 21)
        v = zeros(Float64, 21)

        g₀ = gamma(Complex(1/2, a))
        ϕ₂ = imag(g₀)

        gref = gamma(Complex(1/2, a))
        gr₀, gi₀ = real(gref), imag(gref)
        den = gr₀^2 + gi₀^2

        for k in 2:2:40
            m = k ÷ 2
            g = gamma(Complex(k + 0.5, a))
            gr, gi = real(g), imag(g)
            u[m] = (gr * gr₀ + gi * gi₀) / den
            v[m] = (gr₀ * gi - gr * gi₀) / den
        end

        x² = x^2
        sv₁ = v[1] / (2x²)
        su₂ = u[1] / (2x²)
        fac = 1.0
        for k in 3:2:19
            fac *= -k * (k - 1)
            denom = fac * 2^k * x^(2k)
            sv₁ += v[k] / denom
            su₂ += u[k] / denom
        end

        sv₂ = 0.0
        su₁ = 0.0
        fac = 1.0
        for k in 2:2:20
            fac *= -k * (k - 1)
            denom = fac * 2^k * x^(2k)
            sv₂ += v[k] / denom
            su₁ += u[k] / denom
        end

        s₁ = 1 + sv₁ - su₁
        s₂ = -sv₂ - su₂

        ea = exp(pi * a)
        S = sqrt(1 + ea^2)
        f₊ = S + ea
        f₋ = S - ea

        ϕ = x^2 / 4 - a * log(abs(x)) + pi / 4 + ϕ₂ / 2

        if x > 0
            return sqrt(2 * f₋ / abs(x)) * (s₁ * cos(ϕ) - s₂ * sin(ϕ))
        else
            return sqrt(2 * f₊ / abs(x)) * (s₁ * sin(ϕ) + s₂ * cos(ϕ))
        end
    end
end

"""
    dU(a::Float64, x::Float64)::Float64

Compute the derivative of the parabolic cylinder function U(a,x) for real parameters.
"""
function dU(a::Float64, x::Float64)::Float64
    ε = 1e-15

    c = zeros(Float64, 101)
    c[1] = a
    c₀, c₁ = 1.0, a
    for k in 4:2:202
        m = k ÷ 2
        c[m] = a * c₁ + (k - 2) * (k - 3) * c₀ / 4
        c₀, c₁ = c₁, c[m]
    end

    y₁ = a
    term = 1.0
    @inbounds for k in 1:100
        term *= 0.5 * x^2 / (k * (2k + 1))
        Δ = c[k + 1] * term
        y₁ += Δ
        if abs(Δ / y₁) ≤ ε && k > 30
            break
        end
    end
    y₁ *= x

    d = zeros(Float64, 101)
    d[1], d[2] = 1.0, a
    d₁, d₂ = 1.0, a
    for k in 5:2:160
        m = (k + 1) ÷ 2
        d[m] = a * d₂ + 0.25 * (k - 2) * (k - 3) * d₁
        d₁, d₂ = d₂, d[m]
    end

    y₂ = 1.0
    term = 1.0
    @inbounds for k in 1:100
        term *= 0.5 * x^2 / (k * (2k - 1))
        Δ = d[k + 1] * term
        y₂ += Δ
        if abs(Δ / y₂) ≤ ε && k > 30
            break
        end
    end

    if a < 0 && isapprox(a + 0.5, round(a + 0.5))
        θ = pi * (0.25 + a / 2)
        f₁ = gamma(0.25 - a / 2) / (√π * 2^(a / 2 + 0.25))
        f₂ = gamma(0.75 - a / 2) / (√π * 2^(a / 2 - 0.25))
        du = cos(θ) * f₁ * y₁ - sin(θ) * f₂ * y₂
    else
        prefactor = √π / 2^(a / 2 + 0.25)
        g₁ = gamma(0.25 + a / 2)
        g₃ = gamma(0.75 + a / 2)
        du = prefactor * (y₁ / g₃ - √2 * y₂ / g₁)
    end

    return du
end

"""
    dV(a::Float64, x::Float64)::Float64

Compute the derivative of the parabolic cylinder function V(a,x) for real parameters.
"""
function dV(a::Float64, x::Float64)::Float64
    ε = 1e-15

    c = zeros(Float64, 101)
    c[1] = a
    c₀, c₁ = 1.0, a
    for k in 4:2:202
        m = k ÷ 2
        c[m] = a * c₁ + (k - 2) * (k - 3) * c₀ / 4
        c₀, c₁ = c₁, c[m]
    end

    y₁ = a
    term = 1.0
    @inbounds for k in 1:100
        term *= 0.5 * x^2 / (k * (2k + 1))
        Δ = c[k + 1] * term
        y₁ += Δ
        if abs(Δ / y₁) ≤ ε && k > 30
            break
        end
    end
    y₁ *= x

    d = zeros(Float64, 101)
    d[1], d[2] = 1.0, a
    d₁, d₂ = 1.0, a
    for k in 5:2:160
        m = (k + 1) ÷ 2
        d[m] = a * d₂ + 0.25 * (k - 2) * (k - 3) * d₁
        d₁, d₂ = d₂, d[m]
    end

    y₂ = 1.0
    term = 1.0
    @inbounds for k in 1:100
        term *= 0.5 * x^2 / (k * (2k - 1))
        Δ = d[k + 1] * term
        y₂ += Δ
        if abs(Δ / y₂) ≤ ε && k > 30
            break
        end
    end

    if a < 0 && isapprox(a + 0.5, round(a + 0.5))
        θ = pi * (0.25 + a / 2)
        f₁ = gamma(0.25 - a / 2) / (√π * 2^(a / 2 + 0.25))
        f₂ = gamma(0.75 - a / 2) / (√π * 2^(a / 2 - 0.25))
        dv = (sin(θ) * f₁ * y₁ + cos(θ) * f₂ * y₂) / gamma(0.5 - a)
    else
        sinπa = sinpi(a)
        g₀ = gamma(0.5 + a)
        p₀ = g₀ / (√π * 2^(a / 2 + 0.25))
        g₁ = gamma(0.25 + a / 2)
        g₃ = gamma(0.75 + a / 2)
        dv = p₀ * (y₁ * (1 + sinπa) / g₃ + √2 * y₂ * (1 - sinπa) / g₁)
    end

    return dv
end



"""
    dW(a::Float64, x::Float64)::Float64

Compute the derivative of the parabolic cylinder function W with parameters `a` evaluated at `x`.
"""
function dW(a::Float64, x::Float64)::Float64
    ε = 1e-15
    p₀ = 2^(-3/4)

    g₁ = abs(gamma(Complex(1/4, a/2)))
    g₃ = abs(gamma(Complex(3/4, a/2)))
    f₁ = sqrt(g₁ / g₃)
    f₂ = sqrt(2 * g₃ / g₁)

    c = zeros(Float64, 101)
    c[1] = a
    c₀, c₁ = 1.0, a
    for k in 4:2:202
        m = k ÷ 2
        c[m] = a * c₁ - (k - 2) * (k - 3) * c₀ / 4
        c₀, c₁ = c₁, c[m]
    end

    y₁ = a
    term = 1.0
    @inbounds for k in 1:100
        term *= 0.5 * x^2 / (k * (2k + 1))
        Δ = c[k + 1] * term
        y₁ += Δ
        if abs(Δ / y₁) ≤ ε && k > 30
            break
        end
    end
    y₁ *= x

    d = zeros(Float64, 101)
    d[1], d[2] = 1.0, a
    d₁, d₂ = 1.0, a
    for k in 5:2:160
        m = (k + 1) ÷ 2
        d[m] = a * d₂ - 0.25 * (k - 2) * (k - 3) * d₁
        d₁, d₂ = d₂, d[m]
    end

    y₂ = 1.0
    term = 1.0
    @inbounds for k in 1:100
        term *= 0.5 * x^2 / (k * (2k - 1))
        Δ = d[k + 1] * term
        y₂ += Δ
        if abs(Δ / y₂) ≤ ε && k > 30
            break
        end
    end

    return p₀ * (f₁ * y₁ - f₂ * y₂)
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

