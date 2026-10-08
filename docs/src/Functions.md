# Functions

This page describes each family of functions: its definition, the accepted
arguments, the method, and examples. Number types and error estimates are
collected in [Accuracy and number types](@ref), derivative support in
[Automatic differentiation](@ref), and the docstrings in the
[API reference](@ref).

The plotting examples use [Plots.jl](https://github.com/JuliaPlots/Plots.jl)
and [LaTeXStrings.jl](https://github.com/JuliaStrings/LaTeXStrings.jl), and
two of them use
[DomainColoring.jl](https://github.com/eprovst/DomainColoring.jl). These
packages are not dependencies of
FewSpecialFunctions.jl; install them with
`import Pkg; Pkg.add(["Plots", "LaTeXStrings", "DomainColoring"])` to
reproduce the figures.

## Coulomb wave functions

The Coulomb wave functions solve the radial Schrödinger equation for a
charged particle in a Coulomb potential,

```math
\frac{d^2u}{d\rho^2} + \left(1 - \frac{2\eta}{\rho} - \frac{\ell(\ell+1)}{\rho^2}\right)u = 0 .
```

All functions share three arguments: the angular momentum ``\ell``, the
Coulomb (Sommerfeld) parameter ``\eta``, which is positive for a repulsive
and negative for an attractive interaction, and the dimensionless radius
``\rho = kr``. Real and complex values are accepted for all three.

- [`F`](@ref)`(ℓ, η, ρ)` is the regular solution, ``F_\ell = C_\ell(\eta)\,\rho^{\ell+1}e^{-i\rho}\,{}_1F_1(\ell+1-i\eta;2\ell+2;2i\rho)``.
- [`G`](@ref)`(ℓ, η, ρ)` is the irregular solution, ``G_\ell = (H^+_\ell + H^-_\ell)/2``.
- [`H⁺`](@ref) and [`H⁻`](@ref) are the outgoing and incoming solutions,
  defined through Tricomi's function ``U``; for real arguments ``H^\pm_\ell = G_\ell \pm iF_\ell``.
- [`C`](@ref)`(ℓ, η)` is the normalization constant ``C_\ell(\eta) = 2^\ell e^{-\pi\eta/2}\,|\Gamma(\ell+1+i\eta)|/\Gamma(2\ell+2)`` for real arguments.
- [`θ`](@ref)`(ℓ, η, ρ)` is the phase ``\rho - \ell\pi/2 - \eta\log 2\rho + \sigma_\ell(\eta)`` that the solutions approach at large ``\rho``.
- [`η`](@ref)`(a, k)` computes the Coulomb parameter ``1/(ak)`` from the Bohr radius and wave number.

`F`, `G`, `C`, and `θ` return real values for real arguments. The remaining
functions ([`D⁺`](@ref), [`D⁻`](@ref), [`Φ`](@ref), [`M_regularized`](@ref),
[`w`](@ref), [`w_plus`](@ref), [`w_minus`](@ref), [`h_plus`](@ref),
[`h_minus`](@ref), [`g`](@ref), [`Φ_dot`](@ref), [`F_dot`](@ref), [`Ψ`](@ref),
and `I`) are the building blocks of the connection formulas in
[arXiv:1804.10976](https://arxiv.org/abs/1804.10976), which the
implementation follows; their docstrings give the formula each one evaluates.

```jldoctest coulomb
julia> using FewSpecialFunctions

julia> ℓ, η0, ρ = 0, 0.3, 2.0;

julia> H⁺(ℓ, η0, ρ) ≈ G(ℓ, η0, ρ) + im * F(ℓ, η0, ρ)
true

julia> isapprox(F(ℓ, η0, 1.0e-4), C(ℓ, η0) * 1.0e-4; rtol = 1.0e-3)
true
```

The confluent hypergeometric functions are evaluated by
[HypergeometricFunctions.jl](https://github.com/JuliaMath/HypergeometricFunctions.jl), which determines the accuracy; see
[Accuracy and number types](@ref).

```@example
using Plots, FewSpecialFunctions, LaTeXStrings # hide
ENV["GKSwstype"] = "100" # hide

plot_font = "Computer Modern" # hide
default(fontfamily=plot_font,linewidth=2.5, framestyle=:box, label=nothing, grid=true,palette=:tab10) # hide
x = range(0, stop=25, length=1000)
plot(x, real(F.(0.0, 0.3, x)), label=L"F_0(0.3,ρ)")
plot!(x, real(F.(0.0, -0.3, x)), label=L"F_0(-0.3,ρ)")
xlabel!(L"ρ")
title!("Regular Coulomb Wave Functions")
```

The same approach shows the regular Coulomb functions for several values of ``\ell``:

```@example
using Plots, FewSpecialFunctions, LaTeXStrings # hide
ENV["GKSwstype"] = "100" # hide

plot_font = "Computer Modern" # hide
default(fontfamily=plot_font,linewidth=2.5, framestyle=:box, label=nothing, grid=true,palette=:tab10) # hide
x = range(0, stop=25, length=1000)
plot(x, real(F.(1e-5, 5.0, x)), label=L"F_0(5.0,ρ)", linewidth=2)
plot!(x, real(F.(1.0, 5.0, x)), label=L"F_1(5.0,ρ)", linewidth=2)
plot!(x, real(F.(2.0, 5.0, x)), label=L"F_2(5.0,ρ)", linewidth=2)
plot!(x, real(F.(3.0, 5.0, x)), label=L"F_3(5.0,ρ)", linewidth=2)
title!("Regular Coulomb Wave Functions for Different ℓ")
xlabel!(L"ρ")
```

### Complex arguments

A domain coloring of ``F_0(z, z)`` in the complex plane:

```@example Complex_Coulomb
using DomainColoring, FewSpecialFunctions, Plots

domaincolor(z -> F(0,z,z), [-2, 2, 0, 5], grid=true)
```

## Whittaker functions

[`WhittakerM`](@ref) and [`WhittakerW`](@ref) evaluate the standard solutions
of Whittaker's equation:

```math
M_{\kappa,\mu}(z) = e^{-z/2}z^{\mu+1/2}
 {}_1F_1(\mu-\kappa+1/2;1+2\mu;z),\qquad
W_{\kappa,\mu}(z) = e^{-z/2}z^{\mu+1/2}
 U(\mu-\kappa+1/2,1+2\mu,z).
```

Here the three-argument ``U`` is Tricomi's hypergeometric function, distinct
from the package's two-argument parabolic cylinder function `U(a, x)`.
[`dWhittakerM`](@ref) and [`dWhittakerW`](@ref) differentiate with respect to `z`.

Inputs must be finite, and this API requires `z ≠ 0`. Positive real arguments
with real parameters return real values. Complex inputs select the principal
branch; use `complex(z)` explicitly for negative real arguments. The sign
of a zero imaginary part selects the corresponding side of the cut.
`WhittakerM` has parameter poles at negative integer `2μ`, which raise
`DomainError`; `WhittakerW` remains defined there.

```jldoctest whittaker
julia> using FewSpecialFunctions

julia> WhittakerM(0, 0.5, 2) ≈ 2sinh(1)
true

julia> WhittakerW(0, 0.5, 2) ≈ exp(-1)
true

julia> dWhittakerW(0, 0.5, 2) ≈ -exp(-1) / 2
true

julia> WhittakerW(0.2, 0.3, 1 - 2im) ≈ conj(WhittakerW(0.2, 0.3, 1 + 2im))
true
```

The implementation uses the existing Kummer series with additional working
precision, a connection formula for `W` with removable parameter poles
evaluated by limits, and an asymptotic expansion accepted only when its terms
reach working precision. These methods follow [DLMF 13.14](https://dlmf.nist.gov/13.14)
and the series techniques in [Thompson & Barnett (1986)](https://doi.org/10.1016/0021-9991(86)90046-X).
It does not implement the full COULCC continued-fraction algorithm.
The extra-precision series favors accuracy over speed, especially for large
parameters or complex arguments. Its precision changes share the cylinder
fallback lock described below.

ForwardDiff uses analytic derivatives in `z` and finite differences for real
`κ` and `μ`; see [Automatic differentiation](@ref).
Tests include independent [mpmath](https://mpmath.org) values; the implementation was also compared
with the direct evaluation path in
[WhittakerCoulomb.jl](https://github.com/banana-bred/WhittakerCoulomb.jl/tree/bab2a17).
The latter shares our hypergeometric dependency, so that comparison alone
does not establish numerical accuracy.

### Real arguments

The following plots hold `μ = 0.3` fixed and vary `κ`. For these parameters,
`M` grows at large positive arguments, while `W` decays. Each panel uses its
own vertical scale. The grid starts above zero because the API requires
a nonzero argument.

```@example WhittakerReal
using Plots, FewSpecialFunctions, LaTeXStrings # hide
ENV["GKSwstype"] = "100" # hide
default(fontfamily = "Computer Modern", linewidth = 2.5, framestyle = :box, grid = true, palette = :tab10) # hide

x = range(0.05, 8, length = 250)
μ = 0.3
pM = plot(xlabel = L"x", ylabel = L"M_{\kappa,0.3}(x)", title = "Whittaker M", legend = :topleft)
pW = plot(xlabel = L"x", ylabel = L"W_{\kappa,0.3}(x)", title = "Whittaker W", legend = :topright)
for (κ, style) in zip((-0.5, 0.0, 0.5), (:solid, :dash, :dot))
    plot!(pM, x, WhittakerM.(κ, μ, x), label = "κ = $κ", linestyle = style)
    plot!(pW, x, WhittakerW.(κ, μ, x), label = "κ = $κ", linestyle = style)
end
plot(pM, pW, layout = (1, 2), size = (800, 350))
```

### Argument derivatives

[`dWhittakerM`](@ref) and [`dWhittakerW`](@ref) give the slopes of the
corresponding solutions. Here both parameters are fixed, and the derivatives
are evaluated directly using the analytic formulas in the package.

```@example WhittakerDerivatives
using Plots, FewSpecialFunctions, LaTeXStrings # hide
ENV["GKSwstype"] = "100" # hide
default(fontfamily = "Computer Modern", linewidth = 2.5, framestyle = :box, grid = true, palette = :tab10) # hide

κ, μ = 0.2, 0.3
x = range(0.1, 4, length = 250)
plot(
    x, dWhittakerM.(κ, μ, x), label = L"M'_{0.2,0.3}(x)",
    xlabel = L"x", ylabel = "derivative", title = "Whittaker argument derivatives",
    legend = :topleft, size = (650, 350)
)
plot!(x, dWhittakerW.(κ, μ, x), label = L"W'_{0.2,0.3}(x)", linestyle = :dash)
```

### Complex arguments

This example follows the vertical line `z = 1 + iy`, which avoids the origin
and the negative-real branch cut. For real `κ` and `μ`, conjugation symmetry
makes the real parts even in `y` and the imaginary parts odd.

```@example WhittakerComplex
using Plots, FewSpecialFunctions, LaTeXStrings # hide
ENV["GKSwstype"] = "100" # hide
default(fontfamily = "Computer Modern", linewidth = 2.5, framestyle = :box, grid = true, palette = :tab10) # hide

κ, μ = 0.2, 0.3
y = range(-6, 6, length = 251)
z = 1 .+ im .* y
m = WhittakerM.(κ, μ, z)
w = WhittakerW.(κ, μ, z)
pM = plot(
    y, real.(m), label = "real part", xlabel = L"y", ylabel = "value",
    title = L"M_{0.2,0.3}(1+iy)", legend = :outerbottom, legend_columns = 2
)
plot!(pM, y, imag.(m), label = "imaginary part", linestyle = :dash)
pW = plot(
    y, real.(w), label = "real part", xlabel = L"y", ylabel = "value",
    title = L"W_{0.2,0.3}(1+iy)", legend = :outerbottom, legend_columns = 2
)
plot!(pW, y, imag.(w), label = "imaginary part", linestyle = :dash)
plot(pM, pW, layout = (1, 2), size = (800, 400))
```

## Marcum Q-function

The Marcum Q-function is a generalized integral involving the modified Bessel function of the first kind. It is widely used in communications and radar signal processing. The implementation in this package is based on the methods described in [arXiv:1311.0681v1](https://arxiv.org/pdf/1311.0681v1), providing accurate results for a wide range of parameters.

- `MarcumQ(μ, a, b)`: Computes the Marcum Q-function for order `μ`, non-centrality parameter `a`, and threshold `b`. The implementation uses series expansions, asymptotic expansions, and recurrence relations for efficiency and accuracy.
- `dQdb(M, a, b)`: Computes the derivative of the Marcum Q-function with respect to `b`.

The code automatically selects the most appropriate algorithm depending on the input parameters.

```@example MarcumQ
using Plots, FewSpecialFunctions, LaTeXStrings  # hide
ENV["GKSwstype"] = "100" # hide

plot_font = "Computer Modern" # hide
default(fontfamily=plot_font,linewidth=2.5, framestyle=:box, label=nothing, grid=true,palette=:tab10) # hide
bs = collect(range(0.0,10,length=100))
M1 = MarcumQ(1, 0.2, bs)
M2 = MarcumQ(1, 1.3, bs)
M3 = MarcumQ(1, 2.5, bs)
M4 = MarcumQ(1, 4.7, bs)

plot(bs, M1, label=L"a=0.2")
plot!(bs, M2, label=L"a=1.3")
plot!(bs, M3, label=L"a=2.5")
plot!(bs, M4, label=L"a=4.7")
plot!(xlabel="b", ylabel=L"Q(1,a,b)", title="Marcum Q-function")
```

### Derivative of the Marcum Q-function

```@example MarcumQ_derivative
using Plots, FewSpecialFunctions, LaTeXStrings  # hide
ENV["GKSwstype"] = "100" # hide

plot_font = "Computer Modern" # hide
default(fontfamily=plot_font,linewidth=2.5, framestyle=:box, label=nothing, grid=true,palette=:tab10) # hide

bs = collect(range(0.0,10,length=100))
M1 = dQdb(1, 0.2, bs)
plot(bs, M1, label=L"a=0.2")
```

## Voigt function

The real Voigt function `voigt(x, y)` is the convolution of a Gaussian and a
Lorentzian profile, with nonnegative width parameter `y`.

```@example Voigt
using Plots, FewSpecialFunctions, LaTeXStrings # hide
ENV["GKSwstype"] = "100" # hide

default(fontfamily="Computer Modern", linewidth=2.5, framestyle=:box, grid=true)
x = range(-6, 6, length=1000)
plot(x, voigt.(x, 0.1), label=L"y=0.1", xlabel=L"x", ylabel=L"K(x,y)",
     title="Voigt function")
plot!(x, voigt.(x, 1.0), label=L"y=1.0")
```

## Parabolic cylinder functions

The parabolic cylinder functions `U(a, x)` and `V(a, x)` solve the parabolic cylinder differential equation. `U` and `W` use asymptotic expansions only when their decreasing terms reach the target precision. Otherwise, they use a convergent series with extra working precision to control cancellation. This fallback favors accuracy over speed. The derivatives `dU`, `dV`, and `dW` differentiate the same expansions as their corresponding functions.

- `U(a, x)`: Computes the parabolic cylinder function of the first kind using a combination of series and asymptotic expansions.
- `V(a, x)`: Computes the second, linearly independent solution using a convergent series with extra working precision.
- `W(a, x)`: Solves the related equation ``W''=(a-x^2/4)W`` using series or asymptotic expansions.

The formulas follow [DLMF Chapter 12](https://dlmf.nist.gov/12), including its [expansions for W](https://dlmf.nist.gov/12.14).

On Julia versions with process-global `BigFloat` precision (including Julia 1.10), the extra-precision fallback is serialized between cylinder calls. Avoid running it concurrently with unrelated `BigFloat` arithmetic, which shares that precision setting.

### The ``D_ν`` convention and scaled values

[`ParabolicCylinderD`](@ref) implements ``D_ν(x)=U(-ν-1/2,x)`` for real
order and argument. [`dParabolicCylinderD`](@ref) gives its argument derivative.

For finite real `a` and `x ≥ 0`, [`U_scaled`](@ref) and [`V_scaled`](@ref)
use the scaling of [Gil, Segura & Temme (2006), equations (7), (11)–(13)](https://ir.cwi.nl/pub/14654/14654D.pdf):

```math
U_{\rm scaled}(a,x)=e^{L(a,x)}U(a,x),\qquad
V_{\rm scaled}(a,x)=e^{-L(a,x)}V(a,x),
```

where, with ``q=x^2/4+a``,

```math
L(a,x)=\begin{cases}
x^2/4,&a=0,\\
\frac{a}{2}(\log|a|-1),&q\le0,\ a\ne0,\\
a\log(x/2+\sqrt{q})+\frac{x}{2}\sqrt{q}-\frac{a}{2},&q>0,\ a\ne0.
\end{cases}
```

[`ParabolicCylinderD_scaled`](@ref) is `U_scaled(-ν-1/2, x)`.
The scaling removes dominant behavior in both order and argument; for general
order it differs from simply multiplying `Dν(x)` by `exp(x²/4)`.
The asymptotic path combines exponential factors algebraically; the series
path scales before conversion to the output type. Thus intermediate underflow
or overflow does not destroy the scaled result. Large-order series can be slow;
this implementation does not port Algorithm 850's full method-selection scheme.

```jldoctest cylinder_scaled
julia> using FewSpecialFunctions

julia> ParabolicCylinderD(1, 2) ≈ 2exp(-1)
true

julia> U(0, 60) == 0 && isinf(V(0, 60))
true

julia> U_scaled(0, 60) ≈ 0.12908600517683427
true

julia> V_scaled(0, 60) ≈ 0.10301719023914642
true
```

All five additions preserve `Float32`, `Float64`, and `BigFloat` and support
broadcasting. ForwardDiff supports derivatives in real order and argument;
for the scaled functions the argument derivative includes the derivative
of the scaling factor; see [Automatic differentiation](@ref).

```@example U
using Plots, FewSpecialFunctions, LaTeXStrings # hide
ENV["GKSwstype"] = "100" # hide

plot_font = "Computer Modern" # hide
default(fontfamily=plot_font,linewidth=2.5, framestyle=:box, label=nothing, grid=true,palette=:tab10) # hide
xs = collect(-2.5:0.05:2.5)
plot(xs, U(0.5, xs), label=L"a=0.5")
plot!(xs, U(2.0, xs), label=L"a=2.0")
plot!(xs, U(3.5, xs), label=L"a=3.5")
plot!(xs, U(5.0, xs), label=L"a=5.0")
plot!(xs, U(8.0, xs), label=L"a=8.0")
plot!(xlabel="x", ylabel=L"U(a,x)", title="Parabolic Cylinder Function U(a,x)")
```

### `V(a, x)`
The function `V(a, x)` is a second, linearly independent solution to the same differential equation satisfied by `U(a, x)`. It uses a convergent series with extra working precision.

```@example Cylinder
using Plots, FewSpecialFunctions, LaTeXStrings  # hide
ENV["GKSwstype"] = "100" # hide

plot_font = "Computer Modern" # hide
default(fontfamily=plot_font,linewidth=2.5, framestyle=:box, label=nothing, grid=true,palette=:tab10) # hide
xs = collect(range(-2.5,2.5,length=100))
V1 = V(0.5, xs)
V2 = V(2.0, xs)
V3 = V(3.5, xs)
V4 = V(5.0, xs)

plot(xs, V1, label=L"a=0.5")
plot!(xs, V2, label=L"a=2.0")
plot!(xs, V3, label=L"a=3.5")
plot!(xs, V4, label=L"a=5.0")
ylims!(-3.0, 3.0)  
plot!(xlabel="x", ylabel=L"V(a,x)", title="Parabolic cylinder function V(a,x)")
```


## Debye functions

The Debye functions are given by

```math
    D_n(\beta,x)= \frac{n}{x^n} \int_0^x \frac{t^n}{(\text{e}^t-1)^\beta} \, \text{d}t
```
For finite `n > 0`, the integral converges when `0 < β < n+1`. At `x = 0`, the limiting value is `0` for `β < 1`, `1` for `β = 1`, and `Inf` for `β > 1`.

The implementation uses transformed adaptive quadrature. `tol` specifies a relative tolerance, floored at eight times machine epsilon; `max_terms` limits the number of quadrature subintervals. Failure to reach that tolerance raises an error.

```@example
using Plots, FewSpecialFunctions, LaTeXStrings # hide
ENV["GKSwstype"] = "100" # hide

plot_font = "Computer Modern" # hide
default(fontfamily=plot_font,linewidth=2.5, framestyle=:box, label=nothing, grid=true,palette=:tab10) # hide
x = range(0,stop=25,length=1000)
plot(x,debye_function(1.0,1.0,x),label=L"D_1(x)")
plot!(x,debye_function(2.0,1.0,x),label=L"D_2(x)")
plot!(x,debye_function(3.0,1.0,x), label=L"D_3(x)")
title!("Debye Functions")
xlabel!(L"x")
savefig("debye.svg"); nothing # hide
```
![](./debye.svg)

## Fermi-Dirac integrals

In solid state physics the Fermi-Dirac integral is given by

```math
    F_j(x) = \int_0^\infty \frac{t^j}{\exp(t-x)+1} \, dt.
```
[`FermiDiracIntegral`](@ref)`(j, x)` approximates this unnormalized integral,
and [`FermiDiracIntegralNorm`](@ref)`(j, x)` divides it by ``\Gamma(j+1)``.
Both require ``j \ge -1/2`` and throw an error for smaller orders.

The method, and with it the accuracy, depends on the order:

| Order ``j`` | Method | Measured relative error |
|:--|:--|:--|
| ``0`` | closed form ``\log(1+e^x)`` | ``1.5\times10^{-14}`` for ``x \ge -5`` |
| ``-1/2,\ 1/2,\ 3/2,\ 5/2`` | Antia's rational approximations, in ``e^x`` for ``x<2`` and in ``1/x^2`` for ``x \ge 2`` | ``2.8\times10^{-12}`` or less |
| any other ``j > -1/2`` | [closed-form expression of Aymerich-Humet, Serra-Mestres, and Millan](https://doi.org/10.1063/1.332276) | ``5.8\times10^{-3}`` to ``1.2\times10^{-2}`` at the tested orders ``1/4, 1, 2, 3, 9/2``; ``2.9\times10^{-2}`` at order ``10`` |

Each figure is the maximum over a grid in ``-30 \le x \le 60`` relative to
high-precision quadrature; the per-order values and the grid are given in
[Accuracy and number types](@ref). The errors are properties of the
approximations, not of the number type: `BigFloat` arguments return a
`BigFloat` with the same error. For ``j = 0`` and ``x \ll 0`` the relative
error grows (``1.0\times10^{-3}`` at ``x=-30``) while the absolute error
stays below machine epsilon, and the result overflows to `Inf` once ``e^x``
does.

```jldoctest fermi_dirac
julia> using FewSpecialFunctions

julia> FermiDiracIntegral(3 / 2, 1.0) ≈ 2.6616826247307124
true

julia> FermiDiracIntegralNorm(1 / 2, 1.0) ≈ FermiDiracIntegral(1 / 2, 1.0) / (sqrt(π) / 2)
true

julia> FermiDiracIntegral(0, 0.0) ≈ log(2)
true
```

The normalized integrals of the four half-integer orders:

```@example
using Plots, FewSpecialFunctions, LaTeXStrings # hide
ENV["GKSwstype"] = "100" # hide

plot_font = "Computer Modern" # hide
default(fontfamily=plot_font,linewidth=2.5, framestyle=:box, label=nothing, grid=true,palette=:tab10) # hide
x = range(0,stop=100,length=10000)
plot(x,FermiDiracIntegralNorm.(-1/2,x),label=L"F_{-1/2}(x)")
plot!(x,FermiDiracIntegralNorm.(1/2,x),label=L"F_{1/2}(x)")
plot!(x,FermiDiracIntegralNorm.(3/2,x),label=L"F_{3/2}(x)")
plot!(x,FermiDiracIntegralNorm.(5/2,x),label=L"F_{5/2}(x)")
xlabel!(L"x")
title!("Fermi-Dirac Integral")
```

## Bose–Einstein integrals

[`BoseEinsteinIntegral`](@ref) evaluates the unnormalized complete integral

```math
    \mathcal{B}_k(\eta) = \int_0^\infty \frac{t^k}{\exp(t-\eta)-1}\,dt,
```

for integer or half-integer orders `k > -1` and real `η ≤ 0`.
[`BoseEinsteinIntegralNorm`](@ref) evaluates its normalized form

```math
    B_k(\eta) = \frac{\mathcal{B}_k(\eta)}{\Gamma(k+1)}
              = \operatorname{Li}_{k+1}(e^\eta).
```

The normalized function also accepts integer or half-integer orders down to
`k = -9/2`. For `k ≤ -1`, it is defined by repeated differentiation,
``dB_k/d\eta = B_{k-1}``; the integral above does not converge at those orders.

Both functions return `0` at `η = -Inf`. At `η = 0`, the normalized value is
``\zeta(k+1)`` for `k > 0`, and the unnormalized value includes the factor
``\Gamma(k+1)``. For supported `k ≤ 0`, the limit as `η → 0⁻` is `Inf`.
Positive `η`, `NaN`, and unsupported orders raise `DomainError`.

```jldoctest bose_einstein
julia> using FewSpecialFunctions

julia> BoseEinsteinIntegralNorm(0.5, -1.0) ≈ 0.4284407345998379
true

julia> BoseEinsteinIntegral(1.5, -1.0) ≈ (3 * sqrt(π) / 4) * 0.3957280103803376
true

julia> BoseEinsteinIntegralNorm(1, 0) ≈ π^2 / 6
true

julia> BoseEinsteinIntegralNorm(-1, -1.0) ≈ inv(expm1(1.0))
true

julia> all(isfinite, BoseEinsteinIntegralNorm.(0.5, [-2.0, -1.0, 0.0]))
true
```

For `Float32` and `Float64`, the implementation uses minimax rational
approximations with convergent-series fallbacks, based on Tables A.5–A.48 of
[Fukushima's 2020 preprint](https://doi.org/10.13140/RG.2.2.21720.65283)
for half-integer orders `-9/2:1:39/2` and integer orders `1:19`.
The approximations target near-double-precision accuracy for `Float64`;
`Float32` results retain single precision. Nonpositive integer orders use
elementary formulas, and higher orders use convergent series. `BigFloat`
evaluation uses convergent series at the working precision.

ForwardDiff differentiates with respect to `η` using the order-lowering
identity. Differentiation with respect to the discrete order `k` raises
`DomainError`.

## Clausen functions

The Clausen functions are implemented for orders 1 through 6, using a combination of series summation and analytic continuation. The code is based on the methods described in [this paper](https://doi.org/10.1007/s10543-023-00944-4).

- `Clausen(n, θ)`: Computes the Clausen function of order `n` at angle `θ`.
- `F_clausen`, `f_n`, `Ci_complex`: Auxiliary functions for advanced use and analytic continuation.

```@example Clausen
using Plots, FewSpecialFunctions, LaTeXStrings # hide
ENV["GKSwstype"] = "100" # hide

plot_font = "Computer Modern" # hide
default(fontfamily=plot_font,linewidth=2.5, framestyle=:box, label=nothing, grid=true,palette=:tab10) # hide

θ = collect(range(0, 2π, length=100))

plot(θ,Clausen.(1, θ), label=L"C_1(θ)")
plot!(θ,Clausen.(2, θ), label=L"C_2(θ)")
plot!(θ,Clausen.(3, θ), label=L"C_3(θ)")
plot!(θ,Clausen.(4, θ), label=L"C_4(θ)")

xlabel!(L"\theta")
title!("Clausen Functions")
```

## Fresnel integrals

The Fresnel integrals are

```math
C(z) = \int_0^z \cos\!\left(\frac{\pi t^2}{2}\right)\,\mathrm{d}t,
\qquad
S(z) = \int_0^z \sin\!\left(\frac{\pi t^2}{2}\right)\,\mathrm{d}t.
```

`fresnel(z)` returns `(C(z), S(z), C(z) + im * S(z))`. The convenience
functions `FresnelC`, `FresnelS`, and `FresnelE` return the corresponding
components. Real and complex arguments are supported; integer arguments are
promoted to floating point.

```jldoctest fresnel
julia> using FewSpecialFunctions

julia> C1, S1, E1 = fresnel(1.0);

julia> round(C1; digits = 10), round(S1; digits = 10)
(0.7798934004, 0.4382591474)

julia> E1 == C1 + im * S1
true

julia> FresnelE(1 + im) ≈ FresnelC(1 + im) + im * FresnelS(1 + im)
true
```

`BigFloat` arguments are supported, with a gap for real arguments of
intermediate magnitude described in [Accuracy and number types](@ref).

### Real axis

```@example FresnelReal
using Plots, FewSpecialFunctions, LaTeXStrings # hide
ENV["GKSwstype"] = "100" # hide

default(fontfamily="Computer Modern", linewidth=2.5, framestyle=:box, grid=true)
x = range(-6, 6, length=1000)
plot(x, FresnelC.(x), label=L"C(x)", xlabel=L"x", ylabel="value",
     title="Fresnel integrals on the real axis")
plot!(x, FresnelS.(x), label=L"S(x)")
```

### Euler spiral

```@example EulerSpiral
using Plots, FewSpecialFunctions, LaTeXStrings # hide
ENV["GKSwstype"] = "100" # hide

default(fontfamily="Computer Modern", linewidth=2.5, framestyle=:box, grid=true)
x = range(-12, 12, length=1500)
plot(FresnelC.(x), FresnelS.(x), label=nothing, xlabel=L"C(x)",
     ylabel=L"S(x)", title="Euler spiral")
```

### Complex plane

```@example FresnelComplex
using DomainColoring, FewSpecialFunctions

domaincolor(FresnelE, [-2, 2, -2, 2], grid=true)
```

The implementation combines a convergent series near zero, stable intermediate
evaluation, and asymptotic expansions for large arguments. It is adapted from [this paper](https://doi.org/10.1007/s11075-023-01654-2)
and the [upstream Fortran source repository](https://github.com/mofrehzaghloul/Fresnel_Integrals).

## Dawson integral

```math
D(x) = e^{-x^2}\int_0^x e^{t^2}\,\mathrm{d}t.
```

[`dawson`](@ref)`(x)` evaluates the Dawson integral for real `x`. It is odd,
is zero at the origin, and approaches `1 / (2x)` for large `|x|`. It is
related to the imaginary error function by
``D(x) = \tfrac{\sqrt{\pi}}{2}e^{-x^2}\operatorname{erfi}(x)``.

```jldoctest dawson
julia> using FewSpecialFunctions

julia> round(dawson(1.0); digits = 10)
0.5380795069

julia> dawson(-0.5) == -dawson(0.5)
true

julia> dawson(big"1.0") isa BigFloat
true
```

```@example Dawson
using Plots, FewSpecialFunctions, LaTeXStrings # hide
ENV["GKSwstype"] = "100" # hide

default(fontfamily="Computer Modern", linewidth=2.5, framestyle=:box, grid=true)
x = range(-8, 8, length=1000)
plot(x, dawson.(x), label=L"D(x)", xlabel=L"x", ylabel="value", title="Dawson integral")
```

The implementation uses the adaptive series and continued fractions described in
[Zaghloul (2023)](https://doi.org/10.1007/s11075-023-01608-8) and is informed by
the [upstream MIT-licensed Fortran implementation](https://github.com/mofrehzaghloul/Dawson).
