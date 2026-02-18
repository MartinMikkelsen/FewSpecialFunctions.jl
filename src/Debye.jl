using SpecialFunctions

export debye_function

"""
    debye_function(n::T, β::T, x::T; tol=T(1e-35), max_terms=2000) where {T <: AbstractFloat}

Compute the generalized Debye function with parameters `n`, `β`, and `x`.
Supports any `AbstractFloat` type (e.g., `Float32`, `Float64`, `BigFloat`).

References:
- [Debye function](https://en.wikipedia.org/wiki/Debye_function)
- [Paper](https://doi.org/10.1007/s10765-007-0256-1)
"""
function debye_function(n::T, β::T, x::T; tol = T(1.0e-35), max_terms = 2000) where {T <: AbstractFloat}
    sum = zero(T)
    c = one(T)

    γn1 = gamma(n + one(T))
    @inbounds for i in 0:max_terms
        ψ = β + i
        P, _ = gamma_inc(n + one(T), ψ * x)
        γlower = P * γn1

        term = c * γlower / ψ^(n + one(T))
        sum += term

        if abs(term) < tol * abs(sum)
            return n * sum / x^n
        end

        c *= (β + i) / (i + 1)
    end

    return n * sum / x^n
end

function debye_function(n::Real, β::Real, x::Real; tol = 1.0e-35, max_terms = 2000)
    T = float(promote_type(typeof(n), typeof(β), typeof(x)))
    return debye_function(T(n), T(β), T(x); tol = T(tol), max_terms = max_terms)
end

function debye_function(n::Real, β::Real, x::AbstractArray{<:Real}; tol = 1.0e-35, max_terms = 2000)
    return [debye_function(n, β, xᵢ; tol = tol, max_terms = max_terms) for xᵢ in x]
end

function debye_function(n::Real, β::AbstractArray{<:Real}, x::Real; tol = 1.0e-35, max_terms = 2000)
    return [debye_function(n, βᵢ, x; tol = tol, max_terms = max_terms) for βᵢ in β]
end

function debye_function(n::AbstractArray{<:Real}, β::Real, x::Real; tol = 1.0e-35, max_terms = 2000)
    return [debye_function(nᵢ, β, x; tol = tol, max_terms = max_terms) for nᵢ in n]
end

function debye_function(β::Real, x::Real; tol = 1.0e-35, max_terms = 2000)
    return debye_function(one(float(promote_type(typeof(β), typeof(x)))), β, x; tol = tol, max_terms = max_terms)
end
