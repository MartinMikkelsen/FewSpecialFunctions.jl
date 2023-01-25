function Debye_function(n,x)
    D = n/(x^n) * quadgk(t -> t^n/(exp^(t)-1),0,x)
    return D    
end
export Debye_function