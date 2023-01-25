#Only for ν larger than -1/2
function Struve(ν,z)
    H_out = (2*(0.5*z)^ν)/(sqrt(π)*gamma(ν+0.5))
    H_int = quadgk(t -> (1-t^2)^(ν-0.5)*sin(t*z),0,1)[1]
    return H_out*H_int
end
#Only for ν larger than -1/2
function Struve_alt(ν,z)
    H_int = (2*(0.5*z)^ν)/(sqrt(π)*gamma(ν+0.5))
    H_out = quadgk(θ -> sin(z*cos(θ))*sin(θ)^(2*ν),0,π/2)[1]
    return H_int*H_out 
end

export Struve, Struve_alt