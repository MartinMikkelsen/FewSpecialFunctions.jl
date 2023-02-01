@doc raw"""
    Debye_function(n,x,min_tol=1e-5)

The Debye function(n,x) given by
    
``D_n(x) = \frac{n}{x^n} \int_0^x \frac{t^n}{e^{t}-1} dx``

Returns the value ``D(n,x)``
"""
function Debye_function(n,x,min_tol=1e-5)
    D = n./(x.^n).*quadgk(t -> (t.^n)/(exp.(t)-1),min_tol,x)[1]
    return D    
end
export Debye_function