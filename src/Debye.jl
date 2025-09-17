using SpecialFunctions

export debye_function

"""
    debye_function(n::Float64, β::Float64, x::Float64; 
                  tol=1e-35, max_terms=2000) -> Float64

Compute the generalized Debye function with parameters `n`, `β`, and `x`.

References:
- [Debye function](https://en.wikipedia.org/wiki/Debye_function)    
- [Paper](https://doi.org/10.1007/s10765-007-0256-1)
"""
function debye_function(n::Float64, β::Float64, x::Float64; tol = 1.0e-35, max_terms = 2000)
    sum = 0.0
    c = 1.0

    γn1 = gamma(n + 1)
    @inbounds for i in 0:max_terms
        ψ = β + i
        P, _ = gamma_inc(n + 1, ψ * x)
        γlower = P * γn1

        term = c * γlower / ψ^(n + 1)
        sum += term

        if abs(term) < tol * abs(sum)
            return n * sum / x^n
        end

        c *= (β + i) / (i + 1)
    end

    return n * sum / x^n
end

function debye_function(n::Real, β::Real, x::Real; tol = 1.0e-35, max_terms = 2000)
    return debye_function(Float64(n), Float64(β), Float64(x); tol = tol, max_terms = max_terms)
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
    return debye_function(1.0, β, x; tol = tol, max_terms = max_terms)
end
