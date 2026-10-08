# Automatic differentiation

FewSpecialFunctions.jl includes a package extension for
[ForwardDiff.jl](https://github.com/JuliaDiff/ForwardDiff.jl). It is loaded
automatically when both packages are loaded;
[ForwardDiff.jl](https://github.com/JuliaDiff/ForwardDiff.jl) itself is not
installed with
[FewSpecialFunctions.jl](https://github.com/MartinMikkelsen/FewSpecialFunctions.jl):

```julia
import Pkg; Pkg.add("ForwardDiff")
```

The extension defines how each function propagates `ForwardDiff.Dual`
numbers. It does so in one of two ways, and the difference matters for
accuracy:

- **Analytic rule.** The derivative is computed from a known identity, such
  as ``C'(x) = \cos(\pi x^2/2)`` for the Fresnel cosine integral. The
  derivative is then as accurate as the functions in the identity.
- **Finite difference.** Where no identity is implemented, the derivative is
  approximated by a central difference of the function itself. The result is
  an approximation with fewer correct digits than the function value, and it
  has the limitations described [below](@ref fd-limits).

Other automatic differentiation packages are not supported.

## Supported arguments

| Function | Analytic rule | Finite difference | Not differentiable |
|:--|:--|:--|:--|
| `FresnelC`, `FresnelS`, `FresnelE`, `fresnel` | real `z` | | |
| `dawson` | `x` | | |
| `Clausen(n, θ)` | `θ` for `n ≥ 2` | `θ` for `n = 1` | order `n` |
| `BoseEinsteinIntegral`, `BoseEinsteinIntegralNorm` | `η` | | order `k` (throws `DomainError`) |
| `FermiDiracIntegral`, `FermiDiracIntegralNorm` | | `j`, `x` | |
| `MarcumQ(M, a, b)` | `b`, only for `M ≥ 1` and `a ≠ 0` (see [below](@ref marcum-derivative)) | `M`, `a` | |
| `dQdb(M, a, b)` | | `M`, `a`, `b`, for `M ≥ 1` and `a ≠ 0` | |
| `debye_function(n, β, x)` | | `n`, `β`, `x` | |
| `voigt(x, y)` | | `x`, `y` | |
| `U`, `V`, `W`, `dU`, `dV`, `dW` | `x` | `a` | |
| `ParabolicCylinderD`, `dParabolicCylinderD` | `x` | `ν` | |
| `U_scaled`, `V_scaled`, `ParabolicCylinderD_scaled` | `x` | `a`, `ν` | |
| `WhittakerM`, `WhittakerW`, `dWhittakerM`, `dWhittakerW` | `z` | `κ`, `μ` | |
| `η` | all arguments | | |
| `θ(ℓ, η, ρ)` | `ρ` | `ℓ`, `η` | |
| `C`, `D⁺`, `D⁻` | | `ℓ`, `η` | |
| `F`, `G`, `H⁺`, `H⁻`, `F_imag`, `Φ` | | `ℓ`, `η`, `ρ` | |
| `M_regularized` | | all arguments | |
| `w(ℓ, η)` | `η` for half-integer `ℓ` | `η` for `ℓ::Integer` | `ℓ` |
| `w_plus`, `w_minus`, `h_plus`, `h_minus`, `g` | | | all arguments (throw `MethodError`) |

Derivatives are taken with respect to real arguments. The Coulomb functions
`Φ_dot`, `F_dot`, `Ψ`, and `I` are themselves finite differences in `ℓ`;
ForwardDiff can differentiate them, but the result is a finite difference of
a finite difference and loses further digits.

## Examples

The derivative of `FresnelC` uses the analytic rule:

```jldoctest ad
julia> using FewSpecialFunctions, ForwardDiff

julia> x0 = 0.7;

julia> ForwardDiff.derivative(FresnelC, x0) ≈ cos((π / 2) * x0^2)
true

julia> ForwardDiff.derivative(FresnelS, x0) ≈ sin((π / 2) * x0^2)
true
```

For `MarcumQ(M, a, b)`, the derivative in `b` uses [`dQdb`](@ref), while the
derivatives in `a` and `M` are finite differences:

```jldoctest ad
julia> M, a, b = 2.0, 1.5, 3.0;

julia> ForwardDiff.derivative(b -> MarcumQ(M, a, b), b) ≈ dQdb(M, a, b)
true

julia> g = ForwardDiff.gradient(v -> MarcumQ(2.0, v[1], v[2]), [1.5, 3.0]);

julia> g[2] ≈ dQdb(2.0, 1.5, 3.0)
true
```

### [Domain of the Marcum Q derivative](@id marcum-derivative)

The derivative in `b` has a narrower domain than the function. `MarcumQ`
accepts `M ≥ 0.5` and `a ≥ 0`, but `dQdb` requires `M ≥ 1` and `a ≠ 0` and
throws an `AssertionError` otherwise. Differentiating `MarcumQ` in `b` for
`0.5 ≤ M < 1` or at `a = 0` therefore fails, including inside a gradient
that contains `b`. Non-integer orders `M ≥ 1` are accepted. The
finite-difference derivatives in `M` and `a` do not have this restriction,
apart from the boundary behavior described [below](@ref fd-limits).

```jldoctest ad
julia> MarcumQ(0.5, 1.0, 2.0) isa Float64
true

julia> try
           ForwardDiff.derivative(b -> MarcumQ(0.5, 1.0, b), 2.0)
       catch err
           err isa AssertionError
       end
true

julia> ForwardDiff.derivative(a -> MarcumQ(0.5, a, 2.0), 1.0) isa Float64
true
```

### Other examples

For the parabolic cylinder function, the first derivative with respect to `x` is
[`dU`](@ref), and the second follows from the differential equation
``U''(a,x) = (x^2/4 + a)\,U(a,x)``:

```jldoctest ad
julia> a, x0 = 0.5, 1.2;

julia> ForwardDiff.derivative(x -> U(a, x), x0) ≈ dU(a, x0)
true

julia> ForwardDiff.derivative(x -> dU(a, x), x0) ≈ (x0^2 / 4 + a) * U(a, x0)
true
```

The Bose–Einstein integrals are differentiable in `η` through the identity
``dB_k/d\eta = B_{k-1}``. Their order is restricted to integers and
half-integers, so a derivative in `k` is rejected:

```jldoctest ad
julia> ForwardDiff.derivative(η -> BoseEinsteinIntegralNorm(1.5, η), -1.0) ≈ BoseEinsteinIntegralNorm(0.5, -1.0)
true

julia> try
           ForwardDiff.derivative(k -> BoseEinsteinIntegralNorm(k, -1.0), 0.5)
       catch err
           err isa DomainError
       end
true
```

## [Limits of the finite-difference derivatives](@id fd-limits)

A finite-difference derivative of ``f`` at ``x`` is computed as

```math
f'(x) \approx \frac{f(x+h) - f(x-h)}{2h},\qquad
h = \sqrt[3]{\epsilon}\,(|x| + 1),
```

where ``\epsilon`` is the machine epsilon of the argument type; for `Float64`,
``h \approx 6\times10^{-6}\,(|x|+1)``.

**Accuracy.** The truncation error is of order ``h^2 f'''`` and the rounding
error of order ``\delta/h``, where ``\delta`` is the absolute error of the
function values. When the function is accurate to machine precision, the
derivative has a relative error of order `1e-11` in `Float64`. Measured
values are `2.1e-11` for `MarcumQ(1, a, 3)` in `a` at `a = 1.5`, `5.2e-11`
for `Clausen(1, θ)` in `θ` at `θ = 1`, and an absolute error of `4.5e-11` in
the Coulomb Wronskian ``F'G - FG' = 1`` at `ℓ = 0`, `η = 0.3`, `ρ = 2`.

When the function is itself an approximation, its error is amplified by
``1/h``. The Fermi–Dirac integral of order `3/2` has a relative error near
`5e-13`. Its `x`-derivative has a measured relative error of `1.7e-12` at
`x = 1`, and of `3.1e-8` at `x = 2`, where the implementation switches
between two rational approximations. With `Float32` arguments the error at
`x = 1` is `3.1e-6`. For Fermi–Dirac orders evaluated with the closed-form
approximation, the derivative is that of the approximation, whose own
relative error is of order `1e-2`. These figures are reproduced by
`docs/accuracy_checks.jl`; see [Reproducing the measurements](@ref).

**Domain boundaries.** The function is evaluated at ``x \pm h``. If either
point lies outside the domain, the function throws and no derivative is
returned. This happens within ``h`` of a boundary, for example for the
derivative of `MarcumQ` in `a` at `a = 0` or in `M` at `M = 0.5`, of
`debye_function` in `x` at `x = 0`, of `voigt` in `y` at `y = 0`, and of
`FermiDiracIntegral` in `j` at `j = -1/2`.

```jldoctest ad
julia> try
           ForwardDiff.derivative(y -> voigt(0.7, y), 0.0)
       catch err
           err isa DomainError
       end
true
```

**Singularities.** A central difference that straddles a singularity or a
discontinuity returns a finite but meaningless number. `Clausen(1, θ)` near
multiples of ``2\pi`` is an example.

**Higher derivatives.** Nested `ForwardDiff` calls work. Where the first
derivative uses an analytic rule that is itself differentiable analytically
(for example `U` and `dU`, the Fresnel integrals, or `dawson`), the second derivative keeps
that accuracy. Where both levels are finite differences, each level loses
further digits.

**Number types.** The step is chosen from the argument type, so `Float32`
and `BigFloat` arguments work, with the accuracy scaling described above.

If you need a derivative with full accuracy and only a finite-difference rule
is available, use an analytic identity for the function where one exists:
for example ``d\mathcal{F}_j/dx = j\,\mathcal{F}_{j-1}`` for the unnormalized
Fermi–Dirac integral, or the Coulomb wave equation for second derivatives in
``\rho``.
