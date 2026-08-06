
## References

- [Clausen functions: Quadrature processes for efficient calculation of the Clausen functions](https://doi.org/10.1007/s10543-023-00944-4)
- [Coulomb wave functions: Connection formulas between Coulomb wave functions](https://arxiv.org/abs/1804.10976)
- [Debye functions: Calculation of Integer and Noninteger n-Dimensional Debye Functions Using Binomial Coefficients and Incomplete Gamma Functions](https://doi.org/10.1007/s10765-007-0256-1)
- [Fermi-Dirac integrals: Notes on Fermi-Dirac Integrals](https://arxiv.org/abs/0811.0116)
- [Fresnel integrals: Calculation of Fresnel integrals of real and complex arguments up to 28 significant digits](https://doi.org/10.1007/s11075-023-01654-2)
- [Dawson integral: Numerical calculation of Dawson's integral and its imaginary error functions of complex arguments for arbitrary value of the phase angle](https://doi.org/10.1007/s11075-023-01608-8)
- [Marcum Q-function: https://arxiv.org/pdf/1311.0681v1](https://arxiv.org/pdf/1311.0681v1)
- [Parabolic cylinder functions](https://www.scirp.org/reference/referencespapers?referenceid=1112080)

## API

```@docs
# Coulomb wave functions
η
C
θ
F
D⁺
D⁻
H⁺
H⁻
F_imag
G
M_regularized
Φ
w
w_plus
w_minus
h_plus
h_minus
g
Φ_dot
F_dot
Ψ
I

# Debye functions
debye_function

# Fresnel integrals
fresnel
FresnelC
FresnelS
FresnelE

# Clausen functions
Clausen
Ci_complex
f_n
F_clausen

# Fermi-Dirac integrals
FermiDiracIntegral
FermiDiracIntegralNorm

# Marcum Q-function
MarcumQ
dQdb

# Voigt function
voigt

# Parabolic cylinder functions
U
V
W
dU
dV
dW
```

## Dawson integral

`dawson(x::Real)` evaluates the real Dawson integral and preserves `Float32`,
`Float64`, and `BigFloat` inputs. Integer and other real inputs are promoted to
floating point. The implementation is related to the imaginary error function
`erfi(x)`, but this documentation expresses that relation in prose rather than
through a cross-reference.
