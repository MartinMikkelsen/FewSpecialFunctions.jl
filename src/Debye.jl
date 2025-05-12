using SpecialFunctions

export debye_function

"""
    debye_function(n::Float64, β::Float64, x::Float64; 
                  tol=1e-35, max_terms=2000) -> Float64

Compute the generalized Debye function with parameters `n`, `β`, and `x`.

References:
- [Debye function](https://en.wikipedia.org/wiki/Debye_function)    
- [Calculation of Integer and Noninteger n-Dimensional Debye Functions Using Binomial Coefficients and Incomplete Gamma Functions](10.1007/s10765-007-0256-1)
"""
function debye_function(n::Float64, β::Float64, x::Float64; tol=1e-35, max_terms=2000)
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

