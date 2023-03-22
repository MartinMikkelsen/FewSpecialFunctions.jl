using SpecialFunctions

export FermiDiracIntegral

@doc raw"""
    FermiDiracIntegral(j,x)

The Fermi-Dirac integral
    
```math
    F_j(x) = 1/(\Gamma(j+1)) \int_0^\infty \frac{t^j}{\exp(t-x)+1} \, dt
```

Returns the value ``F_j(x)``

Resources:

https://de.wikipedia.org/wiki/Fermi-Dirac-Integral
"""
function FermiDiracIntegral(j,x)
    if j==0
        return log(1+exp(x))
    end
    if j==-1
        return exp(x)/(1+exp(x))
    end
    if j==1/2 && x<1.3
        return  1/(exp(-x)+0.27)
    end
    if j==1/2 && x>=1.3 
        return 4/(3*sqrt(π))*(x^2+π^2/6)^(3/4)
    end 
end