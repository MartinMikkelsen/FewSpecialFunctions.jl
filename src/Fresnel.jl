using SpecialFunctions
export fresnel, S, C, E


"""
    fresnel(z::Number)

Compute the Fresnel integrals S(z) and C(z), and the auxiliary E(z),
using the complex error function:

    zp = (√π/2)*(1 - i)*z
    zm = (√π/2)*(1 + i)*z
    ep = erf(zp)
    em = erf(zm)

Then
    resp = (1 + i)/2 * em
    resn = (1 - i)/2 * ep

and
    S = (resp - resn)/(2i),
    C = (resp + resn)/2,
    E = (1 + i)/2 - resp.

Returns a named tuple `(S=S, C=C, E=E)`.
"""
function fresnel(z::Number)
    zc = complex(z)
    zp = (sqrt(π)/2)*(1 - im)*zc
    zm = (sqrt(π)/2)*(1 + im)*zc
    ep = erf(zp)
    em = erf(zm)
    resp = (1 + im)/2 * em
    resn = (1 - im)/2 * ep
    Sval = (resp - resn)/(2*im)
    Cval = (resp + resn)/2
    Eval = (1 + im)/2 - resp
    return Sval, Cval, Eval
end

"""S(x) -> Fresnel sine integral."""
FresnelS(z::Number) = fresnel(z)[1]
"""C(x) -> Fresnel cosine integral."""
FresnelC(z::Number) = fresnel(z)[2]
"""E(z) -> exponential auxiliary function."""
FresnelE(z::Number) = fresnel(z)[3]

