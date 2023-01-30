"""
    bar(x[, y])

Compute the Bar index between `x` and `y`.

If `y` is unspecified, compute the Bar index between all pairs of columns of `x`.

# Examples
```julia-repl
julia> bar([1, 2], [1, 2])
1
```
"""
function Clausen(x, min_tol=1e-15)
    return real((-quadgk(t -> log(2*Complex(sin(t/2))),min_tol,x)[1]))
end

export Clausen