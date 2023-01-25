using SpecialFunctions 
#using the following representation to test the implementation
function S_integral_pi(x)
    S = quadgk(t -> sin(π/2*t^2),0,x)[1]
    return S    
end

function C_integral_pi(x)
    C = quadgk(t -> cos(π/2*t^2),0,x)[1]
    return C    
end

function S_integral(x)
    S = quadgk(t -> sin(t^2),0,x)[1]
    return S    
end

function C_integral(x)
    C = quadgk(t -> cos(t^2),0,x)[1]
    return C    
end

function S_err(x)
    S = sqrt(π/2)*(1+1im)/4*(erf((1+1im)/(sqrt(2))*x)-1im*erf((1-1im)*x/(sqrt(2))))
    return S    
end

function C_err(x)
    C = sqrt(π/2)*(1-1im)/4*(erf((1+1im)/(sqrt(2))*x)+1im*erf((1-1im)/(sqrt(2))*x))
    return C
end

export S_integral_pi, C_integral_pi ,S_integral, C_integral , S_err, C_err

