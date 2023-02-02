using SpecialFunctions 
"""
    Fresnel_S_integral_pi(x)

The Fresnel function S(z) using the definition in Handbook of Mathematical Functions: Abramowitz and Stegun, where
```math
    S(z) = \int_0^x \cos(\pi t^2/2) dt
```
"""
function Fresnel_S_integral_pi(x)
    S = quadgk(t -> sin(π/2*t^2),0,x)[1]
    return S    
end
"""
    Fresnel_C_integral_pi(x)

The Fresnel function C(z) using the definition in Handbook of Mathematical Functions: Abramowitz and Stegun, where

```math
    C(z) = \int_0^x \sin(pi t^2/2) dt
```
"""
function Fresnel_C_integral_pi(x)
    C = quadgk(t -> cos(π/2*t^2),0,x)[1]
    return C    
end
"""
    Fresnel_S_integral(x)

The Fresnel function S(z) using the definition https://en.wikipedia.org/wiki/Fresnel_integral
```math
    S(z) = \int_0^x \sin(t^2) dt
```
"""
function Fresnel_S_integral(x)
    S = quadgk(t -> sin(t^2),0,x)[1]
    return S    
end
"""
    Fresnel_C_integral(x)

The Fresnel function C(z) using the definition https://en.wikipedia.org/wiki/Fresnel_integral
```math
    C(z) = \int_0^x \cos(t^2) dt
```
"""
function Fresnel_C_integral(x)
    C = quadgk(t -> cos(t^2),0,x)[1]
    return C    
end
"""
    Fresnel_S_err(x)

The Fresnel function S(z) using the definition https://en.wikipedia.org/wiki/Fresnel_integral and the error function.
```math
S(z) = \sqrt{\frac{\pi}{2}} \frac{1+i}{4} \left[ \text{erf} \left(\frac{1+i}{\sqrt{2}}z \right) - i \text{erf} \left(\frac{1-i}{\sqrt{2}}z \right]
```
"""
function Fresnel_S_err(x)
    S = sqrt(π/2).*(1+1im)/4*(erf.((1+1im)/(sqrt(2))*x)-1im*erf.((1-1im)*x/(sqrt(2))))
    return real(S)    
end
"""
    Fresnel_C_err(x)

The Fresnel function C(z) using the definition https://en.wikipedia.org/wiki/Fresnel_integral and the error function.
```math
C(z) = \sqrt{\frac{\pi}{2}} \frac{1-i}{4}\left[ \text{erf}\left(\frac{1+i}{\sqrt{2}}z \right) + i \text{erf}\left(\frac{1-i}{\sqrt{2}}z\right]
```
"""
function Fresnel_C_err(x)
    C = sqrt(π/2)*(1-1im)/4*(erf((1+1im)/(sqrt(2))*x)+1im*erf((1-1im)/(sqrt(2))*x))
    return real(C)
end

export Fresnel_S_integral_pi, Fresnel_C_integral_pi , Fresnel_S_integral, Fresnel_C_integral , Fresnel_S_err, Fresnel_C_err

