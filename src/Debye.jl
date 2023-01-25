function Debye_function(n,x,min_tol=1e-5)
    D = n./(x.^n).*quadgk(t -> (t.^n)/(exp.(t)-1),min_tol,x)[1]
    return D    
end
export Debye_function