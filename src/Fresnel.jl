using SpecialFunctions 
#using the following representation to test the implementation
"""
    Fresnel_S_integral_pi(x)

The Fresnel function S(z) using the definition in Handbook of Mathematical Functions: Abramowitz and Stegun, where

```math
    S(z) = ∫_0^x cos(pi/2t^2) dt
```
"""
function Fresnel_S_integral_pi(x)
    S = quadgk(t -> sin(π/2*t^2),0,x)[1]
    return S    
end
"""
The Fresnel function C(z) using the definition in Handbook of Mathematical Functions: Abramowitz and Stegun, where

```math
    C(z) = ∫_0^x sin(π/2t^2) dt
```
"""
function Fresnel_C_integral_pi(x)
    C = quadgk(t -> cos(π/2*t^2),0,x)[1]
    return C    
end
"""
The Fresnel function S(z) using the definition https://en.wikipedia.org/wiki/Fresnel_integral
```math
    S(z) = ∫_0^x sin(t^2) dt
```
"""
function Fresnel_S_integral(x)
    S = quadgk(t -> sin(t^2),0,x)[1]
    return S    
end
"""
The Fresnel function C(z) using the definition https://en.wikipedia.org/wiki/Fresnel_integral
```math
    C(z) = ∫_0^x cos(t^2) dt
```
"""
function Fresnel_C_integral(x)
    C = quadgk(t -> cos(t^2),0,x)[1]
    return C    
end
"""
The Fresnel function S(z) using the definition https://en.wikipedia.org/wiki/Fresnel_integral and the error function.
"""
function Fresnel_S_err(x)
    S = sqrt(π/2).*(1+1im)/4*(erf.((1+1im)/(sqrt(2))*x)-1im*erf.((1-1im)*x/(sqrt(2))))
    return real(S)    
end
"""
The Fresnel function C(z) using the definition https://en.wikipedia.org/wiki/Fresnel_integral and the error function.
"""
function Fresnel_C_err(x)
    C = sqrt(π/2)*(1-1im)/4*(erf((1+1im)/(sqrt(2))*x)+1im*erf((1-1im)/(sqrt(2))*x))
    return real(C)
end

export Fresnel_S_integral_pi, Fresnel_C_integral_pi , Fresnel_S_integral, Fresnel_C_integral , Fresnel_S_err, Fresnel_C_err

