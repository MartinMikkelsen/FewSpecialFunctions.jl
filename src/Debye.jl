@doc raw"""
    Debye_function(n,x,min_tol=1e-15)

The Debye function(n,x) given by
    
```math
    D_n(x) = \frac{n}{x^n} \int_0^x \frac{t^n}{e^{t}-1} dt
```

Returns the value ``D(n,x)``
"""
function Debye_function(n,x,min_tol=1e-5)
    D = n./(x.^n).*quadgk(t -> (t.^n)/(expm1.(t)),min_tol,x)[1]
    return D    
end
export Debye_function
