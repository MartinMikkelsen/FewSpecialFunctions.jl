@doc raw"""
    Struve(ν,z,min_tol=1e-15)

Returns the Struve function given by
```math
    \mathbf{H}_\nu(z) = \frac{2(z/2)^\nu}{\sqrt{\pi}\Gamma(\nu+1/2)} \int_0^1 (1-t)^{{\nu-1/2}}\sin(zt) \, \text{d}t
```
"""
function Struve(ν::Number,z::Number,min_tol=1e-15)
    if ν >= 1
        H_out = (2*(0.5*z).^ν)/(sqrt(π)*gamma(ν+0.5))
        H_int = quadgk(t -> (1-t.^2)^(ν.-0.5)*sin.(t.*z),min_tol,1)[1]
        return H_out*H_int
    end
    if ν == 0
        H_int = (2*(0.5*z).^ν)/(sqrt(π)*gamma(ν+0.5))
        H_out = quadgk(θ -> sin(z*cos(θ))*sin.(θ)^(2*ν),0,π/2)[1]
        return H_int*H_out 
    end
end

export Struve