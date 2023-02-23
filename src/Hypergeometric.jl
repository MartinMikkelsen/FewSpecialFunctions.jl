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
    if real(b)>real(a)>0
        return gamma(b)/(gamma(a)*gamma(b-a))*quadgk(t -> exp(z*t)*t^(a-1)*(1-t)^(-a+b-1),0,1)[1]
    end
    if real(a)>0 
        return 1/(gamma(a))*quadgk(u -> exp(-u/(1-u))*(u/(1-u))^(a-1)*hypergeometric_0F1(b,z*u/(1-u))*u/(1-u)^2,0,1)[1]
    end
end
@doc raw"""
    confluent_hypergeometric_U(a,b,z)
    
Returns the Kummer confluent hypergeometric function 

```math
    U(a,b,z) = \frac{\Gamma(b-1)}{\Gamma(a)}z^{1-b} {}_1 F_1(a-b+1,2-b,z)+\frac{\Gamma(1-b)}{\Gamma(a-b+1)} {}_1F_1(a,b,z)
```
"""
function confluent_hypergeometric_U(a,b,z)
    if real(z)>0 && real(a) > 0 
        return 1/gamma(a)*quadgk(u -> exp(-z*u/(1-u))*(u/(1-u))^(a-1)*(u/(1-u)+1)^(-a+b-1)*u/(1-u).^2,0,1)[1]
    end
end

export hypergeometric_0F1, confluent_hypergeometric_1F1, confluent_hypergeometric_U
