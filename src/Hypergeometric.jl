using SpecialFunctions

@doc raw"""
    hypergeometric_0F1(b,z)

Returns the confluent hypergeometric function given by

```math
    {}_0 F_1(a,b) = \sum_{k=0}^\infty \frac{z^k}{(b)_k k!}
```
for the parameters ``a`` and ``b``
"""
function hypergeometric_0F1(b,z)
    if b == 1
        return 1/π .* quadgk(t -> exp(2*sqrt(z)*cos(t)),0,π)[1]
    end
    if real(b)>=0.5
        return (2*gamma(b))/(sqrt(π)*gamma(b-0.5)).*quadgk(t -> (1-t^2)^(b-1.5)*cosh(2*sqrt(z)*t),0,1)[1]
    end 
end
@doc raw"""
    confluent_hypergeometric_1F1(a,b,z)
    
Returns the Kummer confluent hypergeometric function 

```math
    {}_1 F_1(a,b,z) = \sum_{k=0}^{\infty} \frac{(a)_k z^k}{(b)_k  k!}
```
"""
function confluent_hypergeometric_1F1(a,b,z)
    result = 1.0
    term = 1.0
    for n = 1:100
        term *= (a+n-1) / ((b+n-1)*n) * z
        result += term
        if abs(term) < 1e-12 # check if the term is negligible
            break
        end
    end
    return result
end

@doc raw"""
    confluent_hypergeometric_U(a,b,z)
    
Returns the Kummer confluent hypergeometric function 

```math
    U(a,b,z) = \frac{\Gamma(b-1)}{\Gamma(a)}z^{1-b} {}_1 F_1(a-b+1,2-b,z)+\frac{\Gamma(1-b)}{\Gamma(a-b+1)} {}_1F_1(a,b,z)
```
"""
function confluent_hypergeometric_U(a::Float64, b::Float64, z::Float64)
    f1 = F(a, b, z)
    f2 = F(1+a-2, 2-b, z)
    return gamma(1-b)/(gamma(1+a-b))*f1 + gamma(b-1)/(gamma(a))*f2
end


export hypergeometric_0F1, confluent_hypergeometric_1F1, confluent_hypergeometric_U
