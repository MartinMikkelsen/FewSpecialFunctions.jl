"""
Function which returns the value of the Clausen funtion given by

```math
    C_l(\phi) = - \int_0^\phi \log(|2\sin(x/2)|) dx
```
See https://en.wikipedia.org/wiki/Clausen_function
"""
function Clausen(x, min_tol=1e-15)
    return real((-quadgk(t -> log(2*Complex(sin(t/2))),min_tol,x)[1]))
end

export Clausen