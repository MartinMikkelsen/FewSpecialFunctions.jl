using SpecialFunctions 
#using the following representation to test the implementation
function Fresnel_S_integral_pi(x)
    S = quadgk(t -> sin(π/2*t^2),0,x)[1]
    return S    
end

function Fresnel_C_integral_pi(x)
    C = quadgk(t -> cos(π/2*t^2),0,x)[1]
    return C    
end

function Fresnel_S_integral(x)
    S = quadgk(t -> sin(t^2),0,x)[1]
    return S    
end

function Fresnel_C_integral(x)
    C = quadgk(t -> cos(t^2),0,x)[1]
    return C    
end

function Fresnel_S_err(x)
    S = sqrt(π/2).*(1+1im)/4*(erf.((1+1im)/(sqrt(2))*x)-1im*erf.((1-1im)*x/(sqrt(2))))
    return real(S)    
end

function Fresnel_C_err(x)
    C = sqrt(π/2)*(1-1im)/4*(erf((1+1im)/(sqrt(2))*x)+1im*erf((1-1im)/(sqrt(2))*x))
    return real(C)
end

export Fresnel_S_integral_pi, Fresnel_C_integral_pi , Fresnel_S_integral, Fresnel_C_integral , Fresnel_S_err, Fresnel_C_err

