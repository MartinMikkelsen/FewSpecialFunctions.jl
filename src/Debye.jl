
function Debye_function(n,x)
    D = n/(x^n)*quadgk(t -> (t^n)/(exp(t)-1),1e-5,x)[1]
    return D    
end
export Debye_function