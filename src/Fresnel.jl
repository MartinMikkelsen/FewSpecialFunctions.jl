using SpecialFunctions
export fresnel, FresnelS, FresnelC, FresnelE

"""
fresnel(z::Number)

Compute the Fresnel integrals C(z) and S(z), and the auxiliary E(z),
using the complex error function, matching the convention:

    S(x) = ∫₀ˣ sin(π/2 t^2) dt
    C(x) = ∫₀ˣ cos(π/2 t^2) dt

This matches the MATLAB and NIST convention (engineering, not mathematical).

Returns a tuple (C, S, E), where:
  - C = C(z) = Fresnel cosine integral
  - S = S(z) = Fresnel sine integral
  - E = C(z) + i S(z) (auxiliary value)
"""
function fresnel(z::Number)
    # Correct engineering convention (NIST, MATLAB):
    # C(z) + i S(z) = (1 + i)/2 * erf((sqrt(π)/2) * (1 - i) * z)
    w = (sqrt(π)/2) * (1 - im) * z
    F = (1 + im)/2 * erf(w)
    Cval = real(F)
    Sval = imag(F)
    Eval = F
    return Cval, Sval, Eval
end

"""
    FresnelC(z::Number) -> Number

Computes the Fresnel cosine integral C(z) for the given number `z`.

# Arguments
- `z::Number`: The input value (can be real or complex).

# Returns
- `Number`: The value of the Fresnel cosine integral at `z`.
"""
FresnelC(z::Number) = fresnel(z)[1]
"""
    FresnelS(z::Number) -> Number

Computes the Fresnel sine integral S(z) for a given number `z`.

# Arguments
- `z::Number`: The input value (real or complex) at which to evaluate the Fresnel sine integral.

# Returns
- `Number`: The value of the Fresnel sine integral S(z).
"""
FresnelS(z::Number) = fresnel(z)[2]
"""
    FresnelE(z::Number) -> Number

Computes the Fresnel E integral for the given input `z`.

# Arguments
- `z::Number`: The input value (real or complex) at which to evaluate the Fresnel E integral.

# Returns
- `Number`: The value of the Fresnel E integral at `z`.

"""
FresnelE(z::Number) = fresnel(z)[3]