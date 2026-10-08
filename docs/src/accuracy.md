# Accuracy and number types

Two questions are answered separately on this page: which number types a
function accepts and returns, and how many digits of the result are correct.
They are different properties. A function can return a `BigFloat` whose
digits are correct only to double precision, or only to a relative error of
order `1e-2`, when its algorithm is a fixed approximation.

## Number types

Two aspects of the return value are described separately: its precision and
its structure.

### Precision

Arguments are promoted to a common floating-point type, and the result has
the precision of that type. `Float32`, `Float64`, and `BigFloat` are
supported; integers and rationals are converted with `float`, which gives
`Float64`.

```jldoctest types
julia> using FewSpecialFunctions

julia> typeof(debye_function(2.0f0, 1.0f0, 5.0f0))
Float32

julia> typeof(debye_function(big"2.0", big"1.0", big"5.0"))
BigFloat

julia> typeof(U(0, 1))
Float64

julia> typeof(debye_function(2.0f0, 1.0, 5.0f0))
Float64
```

The exceptions are:

- `dQdb` returns `Float64` for `Float32` arguments and throws a `MethodError`
  for `BigFloat` arguments.
- `h_plus`, `h_minus`, and `g` throw a `MethodError` for `BigFloat`
  arguments, because SpecialFunctions.jl has no complex `BigFloat` digamma
  function.
- `Φ_dot`, `F_dot`, `Ψ`, and `I` take their finite-difference step from the
  type of `ℓ`. With an `Integer` order, which `Ψ` and `I` require for
  integer `ℓ`, the step is the `Float64` one, and `Float32` values of `η`
  and `ρ` give a double-precision result: `Float64` for `F_dot` and
  `ComplexF64` for `Φ_dot`, `Ψ`, and `I`.
- `FresnelC`, `FresnelS`, `FresnelE`, and `fresnel` throw a `MethodError` for
  real `BigFloat` arguments of intermediate magnitude; see below.

### Structure

Most functions return a real number for real arguments. The table lists the
structure of the result for each group of functions; "complex" means a
`Complex` number whose parts have the precision described above.

| Functions | Arguments | Result |
|:--|:--|:--|
| Coulomb `η`, `C`, `θ`, `F`, `G`, `M_regularized`, `w`, `g`, `F_dot` | real | real |
| Coulomb `H⁺`, `H⁻`, `F_imag`, `D⁺`, `D⁻`, `Φ`, `Φ_dot`, `Ψ`, `I`, `w_plus`, `w_minus`, `h_plus`, `h_minus` | real | complex |
| Coulomb functions except `g` | any argument complex | complex |
| Coulomb `g` | real or complex | real |
| Whittaker functions | real with `z > 0` | real |
| Whittaker functions | any argument complex | complex |
| `FresnelC`, `FresnelS` | real | real |
| `FresnelC`, `FresnelS` | complex | complex |
| `FresnelE` | real or complex | complex |
| `fresnel` | real or complex | tuple `(C, S, C + im * S)`; for a real argument the first two entries are real and the third is complex |
| Parabolic cylinder, Debye, Dawson, Clausen, Fermi–Dirac, Bose–Einstein, Marcum Q, Voigt | real | real |

```jldoctest types
julia> fresnel(1.0) isa Tuple{Float64, Float64, ComplexF64}
true

julia> FresnelE(1.0f0) isa ComplexF32
true

julia> H⁺(0, 0.3, 2.0) isa ComplexF64, F(0, 0.3, 2.0) isa Float64
(true, true)
```

The domain of each function (for example `η ≤ 0` for the Bose–Einstein
integrals or `y ≥ 0` for `voigt`) is stated in its docstring in the
[API reference](@ref). Arrays are handled with Julia's dot syntax,
`f.(x)`.

## Accuracy by function

The table states what limits the accuracy of each group of functions and
gives a measured error at the stated test point.

- "Working precision" means that the algorithm iterates or selects its method
  until the result is converged for the precision of the result type, so a
  `BigFloat` argument gives correspondingly more digits.
- "Fixed" means that the method has an error that does not decrease with the
  precision of the type.
- "Finite difference" means that the result is a difference quotient whose
  step, and therefore whose error, depends on the precision of the type.

Every figure is a measurement at the listed test point or over the listed
grid, not a bound over the domain. All figures are relative errors,
``|f_{\rm computed} - f| / |f|``, unless marked as absolute. `BigFloat`
figures are for the default precision of 256 bits (``\epsilon \approx 10^{-77}``).
The measurements are produced by the script described under
[Reproducing the measurements](@ref).

| Functions | Method precision | Test point | `Float64` | `BigFloat` |
|:--|:--|:--|:--|:--|
| Coulomb `F`, `G` | working precision, set by [HypergeometricFunctions.jl](https://github.com/JuliaMath/HypergeometricFunctions.jl) | `ℓ = 0`, `η = 0.3`, `ρ = 2` | `1.2e-15` (`F`), `5.2e-15` (`G`) | `1.6e-76` (`F`) |
| Coulomb `Φ_dot`, `F_dot` | finite difference in `ℓ`, error of order `eps(T)^(2/3)` | `ℓ = 1`, `η = 0.3`, `ρ = 2` | `4.5e-11` (`Φ_dot`), `3.4e-10` (`F_dot`) | `1.3e-51` (`Φ_dot`), `5.3e-51` (`F_dot`) |
| Whittaker `W`, `M` | working precision | `κ = 0`, `μ = 1/2`, `z = 2` | `3.4e-17`, `6.7e-17` | exact to 256 bits |
| `ParabolicCylinderD` | working precision | `ν = 1`, `x = 2` | `3.4e-17` | exact to 256 bits |
| `debye_function` | requested tolerance `tol`, not below `8eps(T)` | `n = 3`, `β = 1`, `x = 2`, default `tol = 1e-35` | `1.8e-16` | `1.1e-62` |
| `FresnelC` | working precision | `x = 1.5` | `2.1e-16` | `2.3e-79` |
| `dawson` | working precision, `4eps(T)` stopping rule | `x = 2.5` | `1.8e-16` | `3.9e-77` |
| `Clausen` | fixed: quadrature rule with `N = 10` or `N = 20` nodes | `n = 2`, `θ = π/2` | `1.2e-16` | `2.5e-20` (`N = 10`), `1.9e-21` (`N = 20`) |
| `BoseEinsteinIntegralNorm` | `Float32`/`Float64`: minimax rational approximations; `BigFloat`: convergent series | `k = 1`, `η = -1` | `6.3e-17` | `1.1e-77` |
| `MarcumQ` | working precision | `M = 1`, `a = 1.5`, `b = 3` (`Q = 0.10`) and `b = 9` (`Q = 7.9e-14`) | `3.3e-17`, `7.4e-16` | `8.2e-68`, `1.1e-76` |
| `voigt` | fixed: Fourier expansion designed for double precision | `x = 0.7`, `y = 0.5` | `1.7e-16` | `9.6e-17` |
| `FermiDiracIntegral` | fixed for all orders except `0`; see the next table | | | |

`Ψ` and `I` are linear combinations of `Φ_dot` values and inherit its error.
With `Float32` arguments, the measured errors at the same test points are
between `4e-9` and `2e-7` for `ParabolicCylinderD`, `debye_function`,
`FresnelC`, `dawson`, `Clausen`, `BoseEinsteinIntegralNorm`, `MarcumQ`, and
`voigt`.

### Fermi–Dirac integrals

The method depends on the order `j` and is the same for every number type.
The table gives the maximum relative error of `FermiDiracIntegral(j, x)` in
`Float64` over a grid of 190 points in ``-30 \le x \le 60`` (spacing 1 for
``x \le -3`` and ``x \ge 9``, spacing 0.1 in between), for each order that
was tested. The reference is the defining integral evaluated by adaptive
quadrature in 512-bit arithmetic.

| Order `j` | Method | Maximum relative error | Attained at |
|:--|:--|:--|:--|
| `-1/2` | rational approximation | `2.8e-12` | `x = 33` |
| `1/2` | rational approximation | `5.4e-13` | `x = 14` |
| `3/2` | rational approximation | `5.1e-13` | `x = -25` |
| `5/2` | rational approximation | `2.5e-13` | `x = 37` |
| `0` | closed form `log(1 + exp(x))` | `1.0e-3` | `x = -30` |
| `1/4` | closed-form approximation | `6.4e-3` | `x = 7` |
| `1` | closed-form approximation | `5.8e-3` | `x = 4.3` |
| `2` | closed-form approximation | `6.0e-3` | `x = 2.6` |
| `3` | closed-form approximation | `7.5e-3` | `x = 9` |
| `9/2` | closed-form approximation | `1.2e-2` | `x = 11` |
| `10` | closed-form approximation | `2.9e-2` | `x = 10` |

The four half-integer orders with rational approximations are the only
orders of that kind; every other nonzero order uses the closed-form
approximation, of which the six listed orders are a sample. For order `0`
the error comes from rounding in `1 + exp(x)` and grows as `x` decreases; it
is `1.5e-14` or less for ``x \ge -5``, and the absolute error stays below
machine epsilon.

With 256-bit `BigFloat` arguments at `x = 1`, the relative error is `4.9e-13`
for order `3/2` and `7.1e-4` for order `1`, the same as in `Float64`.

```jldoctest types
julia> x = FermiDiracIntegral(1, big"1.0");

julia> x isa BigFloat
true

julia> exact = big(π)^2 / 6 + 1 / big(2) - sum((-1)^(n + 1) * exp(-big(n)) / big(n)^2 for n in 1:200);

julia> 1.0e-4 < abs(x - exact) / exact < 1.0e-3
true
```

The reference value in this example uses the identity
``\mathcal{F}_1(x) = \pi^2/6 + x^2/2 + \operatorname{Li}_2(-e^{-x})`` for
``x \ge 0``.

### Other limits of `BigFloat` evaluation

- **Voigt function and Clausen functions.** These return a `BigFloat` with
  roughly double precision (Voigt) or about 20 correct digits (Clausen).
- **Fresnel integrals.** For real `BigFloat` arguments of intermediate
  magnitude, evaluation throws a `MethodError` because the method used in
  that range needs a complex error function that SpecialFunctions.jl does not
  provide for `BigFloat`. Of the arguments 2, 3, 4, 5, 5.5, 6, 10, and 40,
  those from 3 to 10 fail; 1.5, 2, and 40 are evaluated to full precision.

## Reproducing the measurements

The figures on this page, the finite-difference figures in
[Automatic differentiation](@ref), and those in the Fermi–Dirac docstrings
are printed by
[`docs/accuracy_checks.jl`](https://github.com/MartinMikkelsen/FewSpecialFunctions.jl/blob/main/docs/accuracy_checks.jl).
Run it from the repository root:

```
julia --project=docs -e 'using Pkg; Pkg.instantiate()'
julia --project=docs docs/accuracy_checks.jl
```

The script states the reference for each measurement: a closed form, a
convergent series, or adaptive quadrature of the defining integral in
512-bit arithmetic. When a change to the package alters a figure, update the
text to the new output.

## Process-global `BigFloat` precision

Parabolic cylinder and Whittaker functions raise the `BigFloat` working
precision internally when a series needs it. On Julia versions where that
precision is a process-global setting (including Julia 1.10), these calls are
serialized by a lock. Avoid running them concurrently with unrelated
`BigFloat` arithmetic, which shares the setting.
