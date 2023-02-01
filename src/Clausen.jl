@doc raw"""
    Clausen(x, min_tol=1e-15)

Computes the Clausen function  

``Cl_2(\phi) = - \int_0^\phi \log|2\sin(x/2)| dx``

Returns ``Cl_2(\phi)``.
"""
function Clausen(x, min_tol=1e-15)
    return real((-quadgk(t -> log(2*Complex(sin(t/2))),min_tol,x)[1]))
end

export Clausen