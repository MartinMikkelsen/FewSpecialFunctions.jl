@doc raw"""
    Struve(ν::Number,z::Number,min_tol=1e-15)
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