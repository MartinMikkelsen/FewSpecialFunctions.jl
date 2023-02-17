# Functions
## Clausen functions
The Clausen function is given by

```math
Cl_2(\phi)=-\int_0^\phi \log|2\sin(x/2)| \, \text{d}x
```

```@example
using Plots, FewSpecialFunctions, LaTeXStrings

plot_font = "Computer Modern"
default(
    fontfamily=plot_font,
    linewidth=2.5, 
    framestyle=:box, 
    label=nothing, 
    grid=true,
    palette=:tab10,
)


x = range(0,15,1000)
xlabel!(L"ϕ")
title!("Clausen function")
plot(x,Clausen.(x), label=L"Cl_2(ϕ)")
savefig("clausen.svg"); nothing # hide
```
![Clausen function](clausen.svg)

## Debye functions

The Debye functions are given by

```math
    D_n(x)= \frac{n}{x^n} \int_0^x \frac{t^n}{\text{e}^t-1} \, \text{d}x
```
And
```@example
using Plots, FewSpecialFunctions, LaTeXStrings

plot_font = "Computer Modern"
default(
    fontfamily=plot_font,
    linewidth=2.5, 
    framestyle=:box, 
    label=nothing, 
    grid=true,
    palette=:tab10,
)

x = range(0,25,1000)

plot(x,Debye_function.(1,x),label=L"D_1(x)")
plot!(x,Debye_function.(2,x),label=L"D_2(x)")
plot!(x,Debye_function.(3,x), label=L"D_3(x)")
title!("Debye Functions")
xlabel!(L"x")
savefig("debye.svg"); nothing # hide
```
![](./debye.svg)


## Regular Coulomb wave functions

The Coulomb wave equation for a charged particle with arbitrary angular momentum and charge is given by 
```math	
    \nabla^2\psi +\left( k^2-\frac{2\mu}{\hbar^2}V(r)\right)\psi = 0,
```
where ``\mu`` is the reduced mass of the system. The radial wave function ``u(r)`` satisfies the following differential equation
```math
	\frac{\text{d}^2 u_\ell}{\text{d}r^2}+\left( k^2-\frac{\ell(\ell+1)}{r^2}-\frac{2\mu}{\hbar^2}\frac{Ze^2}{r}\right)u_\ell=0,
```
where ``Z`` is the product of the charges. Two independent solutions can be found to this equation -- these are called the regular and irregular Coulomb wave functions denoted ``F_\ell(r)`` and ``G_\ell(r)`` respectively. The regular Coulomb wave function ``F_\ell(r)`` is a real function that vanishes at ``r=0`` and the behaviour of the function is described using a parameter $\eta$ which describes how strongly the Coulomb interaction is
```math
	\eta = \frac{Zmc\alpha }{\hbar k},
```
where ``m`` is the mass of the particle, ``k`` is the wave number and ``\alpha`` is the fine structure constant. The solution to is given by
```math
	F_\ell(\eta,kr) = C_\ell (\eta) (kr)^{\ell+1}\text{e}^{-ikr}  {}_1 F_1(\ell+1-i\eta,2\ell+2,2ikr),
```
where ``{}_1F_1(kr)`` is a confluent hypergeometric function and ``C_\ell(\eta)`` is a normalization constant given by 
```math
	C_\ell(\eta) = \frac{2^\ell \text{e}^{-\pi\eta/2}|\Gamma(\ell+1+i\eta)|}{(2\ell+1)!},
```
where ``\Gamma`` is the gamma function. For numerical purposes, it is useful to use the integral [representation of the regular Coulomb wave function](https://dlmf.nist.gov/33.7)
```math
	F_\ell(\eta,\rho) = \frac{\rho^{\ell+1}2^\ell e^{i\rho-(\pi\eta/2)}}{|\Gamma(\ell+1+i\eta)|} \int_0^1 e^{-2i\rho t}t^{\ell+i\eta}(1-t)^{\ell-i\eta} \, \text{d}t.
```

This implementation need the gamma function from [SpecialFunctions.jl](https://github.com/JuliaMath/SpecialFunctions.jl)
```@example
using Plots, FewSpecialFunctions, LaTeXStrings

plot_font = "Computer Modern"
default(
    fontfamily=plot_font,
    linewidth=2.5, 
    framestyle=:box, 
    label=nothing, 
    grid=true,
    palette=:tab10,
)

x = range(0,25,1000)

plot(x,regular_coulomb.(0,0.3,x), label=L"F_0(0.3,ρ)")
plot!(x,regular_coulomb.(0,-0.3,x), label=L"F_0(-0.3,ρ)")
xlabel!(L"ρ")
title!("Regular Coulomb Wave Functions")
```

Use a similar approach to plot the regular Coulomb functions for different a ``\ell``

```@example
using Plots, FewSpecialFunctions, LaTeXStrings

plot_font = "Computer Modern"
default(
    fontfamily=plot_font,
    linewidth=2.5, 
    framestyle=:box, 
    label=nothing, 
    grid=true,
    palette=:tab10,
)

x = range(0,25,1000)

plot(x,regular_coulomb.(1e-5,5,x), label=L"F_0(5,ρ)")
plot!(x,regular_coulomb.(1,5,x), label=L"F_1(5,ρ)")
plot!(x,regular_coulomb.(2,5,x), label=L"F_2(5,ρ)")
plot!(x,regular_coulomb.(3,5,x), label=L"F_3(5,ρ)")
title!("Regular Coulomb Wave Functions")
xlabel!(L"ρ")
```

## Struve functions
The Struve functions are solutions of the non-homogeneous Bessel's differential equation
```math
    x^2 \frac{\text{d}^2 y}{\text{d}x^2} + x \frac{\text{d}y}{\text{d}x}+(x^2-\alpha^2)y = \frac{4(x/2)^{\alpha+1}}{\sqrt{\pi}\Gamma(\alpha+1/2)}
```
The Struve functions are implemented using the following integral representation
```math
    \mathbf{H}_\nu(z) = \frac{2(z/2)^\nu}{\sqrt{\pi}\Gamma(\nu+1/2)} \int_0^1 (1-t)^{{\nu-1/2}}\sin(zt) \, \text{d}t
```
And
```math
    \mathbf{H}_\nu(z) = \frac{2(z/2)^\nu}{\sqrt{\pi}\Gamma(\nu+1/2)} \int_0^{\pi/2} \sin(z\cos(\theta)) \sin^{2\nu}(\theta) \, \text{d}\theta
```
Here is an example
```@example
using Plots, FewSpecialFunctions, LaTeXStrings

plot_font = "Computer Modern"
default(
    fontfamily=plot_font,
    linewidth=2.5, 
    framestyle=:box, 
    label=nothing, 
    grid=true,
    palette=:tab10,
)

x = range(-5,5,1000)

plot(x,Struve.(0,x),label=L"H_0(x)")
plot!(x,Struve.(1,x),label=L"H_1(x)")
plot!(x,Struve.(2,x),label=L"H_2(x)")
plot!(x,Struve.(3,x),label=L"H_3(x)")
plot!(x,Struve.(4,x),label=L"H_4(x)")
plot!(x,Struve.(5,x),label=L"H_5(x)")
xlabel!(L"x")
title!("Struve Functions")
```
## Fresnel functions

```@example
using Plots, FewSpecialFunctions, LaTeXStrings

plot_font = "Computer Modern"
default(
    fontfamily=plot_font,
    linewidth=2.5, 
    framestyle=:box, 
    label=nothing, 
    grid=true,
    palette=:tab10,
)

x = range(-25,25,5000)

plot(x,Fresnel_C_integral.(x),label=L"C(x)")
plot!(x,Fresnel_C_erf.(x), ls=:dash, lw=1.5, label=L"\tilde{C}(x)")
title!("Fresnel Integral")
xlabel!(L"x")
```
and

```@example
using Plots, FewSpecialFunctions, LaTeXStrings

plot_font = "Computer Modern"
default(
    fontfamily=plot_font,
    linewidth=2.5, 
    framestyle=:box, 
    label=nothing, 
    grid=true,
    palette=:tab10,
)

x = range(-25,25,5000)

plot(x,Fresnel_S_integral.(x),label=L"S(x)")
plot!(x,Fresnel_S_erf.(x), ls=:dash, lw=1.5, label=L"\tilde{S}(x)")
title!("Fresnel Integral")
xlabel!(L"x")
```

### Benchmarks

Using the integral approach

```@example bench 
using FewSpecialFunctions,BenchmarkTools

x = range(0,150,1000)

@benchmark Fresnel_C_erf.($x)
```
Using the error function
```@example bench
using FewSpecialFunctions,BenchmarkTools

@benchmark Fresnel_C_integral.($x)
```