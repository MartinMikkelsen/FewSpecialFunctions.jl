# Functions

## Coulomb wave functions

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
where ``\Gamma`` is the gamma function. For numerical purposes, it is useful to use the [integral representation of the regular Coulomb wave function](https://dlmf.nist.gov/33.7)
```math
	F_\ell(\eta,\rho) = \frac{\rho^{\ell+1}2^\ell e^{i\rho-(\pi\eta/2)}}{|\Gamma(\ell+1+i\eta)|} \int_0^1 e^{-2i\rho t}t^{\ell+i\eta}(1-t)^{\ell-i\eta} \, \text{d}t.
```

This implementation requires the gamma function from [SpecialFunctions.jl](https://github.com/JuliaMath/SpecialFunctions.jl)
```@example
using Plots, FewSpecialFunctions, LaTeXStrings # hide
ENV["GKSwstype"] = "100" # hide

plot_font = "Computer Modern" # hide
default(fontfamily=plot_font,linewidth=2.5, framestyle=:box, label=nothing, grid=true,palette=:tab10) # hide
x = range(0,stop=25,length=1000)
plot(x,regular_Coulomb.(0,0.3,x), label=L"F_0(0.3,ρ)")
plot!(x,regular_Coulomb.(0,-0.3,x), label=L"F_0(-0.3,ρ)")
xlabel!(L"ρ")
title!("Regular Coulomb Wave Functions")
```

Use a similar approach to plot the regular Coulomb functions for different a ``\ell``

```@example
using Plots, FewSpecialFunctions, LaTeXStrings # hide
ENV["GKSwstype"] = "100" # hide

plot_font = "Computer Modern" # hide
default(fontfamily=plot_font,linewidth=2.5, framestyle=:box, label=nothing, grid=true,palette=:tab10) # hide
x = range(0,stop=25,length=1000)
plot(x,regular_Coulomb.(1e-5,5,x), label=L"F_0(5,ρ)")
plot!(x,regular_Coulomb.(1,5,x), label=L"F_1(5,ρ)")
plot!(x,regular_Coulomb.(2,5,x), label=L"F_2(5,ρ)")
plot!(x,regular_Coulomb.(3,5,x), label=L"F_3(5,ρ)")
title!("Regular Coulomb Wave Functions")
xlabel!(L"ρ")
```
## Marcum Q-function
The generalized [Marcum Q-function](https://en.wikipedia.org/wiki/Marcum_Q-function#) is defined as
```math
    Q_\nu (a,b) = \int_b^\infty x^\nu \exp\left(-\frac{x^2+a^2}{2}\right) I_{\nu-1}(ax) \, \text{d}x,
```
where ``b \geq 0``, ``a,\nu>0`` and ``I_{\nu-1}`` is the modified Bessel function of the first kind. 

```@example MarcumQ
using Plots, FewSpecialFunctions, LaTeXStrings  # hide
ENV["GKSwstype"] = "100" # hide

plot_font = "Computer Modern" # hide
default(fontfamily=plot_font,linewidth=2.5, framestyle=:box, label=nothing, grid=true,palette=:tab10) # hide
bs = collect(range(0.0,10,length=100))
M1 = marcum_Q(1, 0.2, bs)
M2 = marcum_Q(1, 1.3, bs)
M3 = marcum_Q(1, 2.5, bs)
M4 = marcum_Q(1, 4.7, bs)

plot(bs, M1, label=L"a=0.2")
plot!(bs, M2, label=L"a=1.3")
plot!(bs, M3, label=L"a=2.5")
plot!(bs, M4, label=L"a=4.7")
plot!(xlabel="b", ylabel=L"Q(1,a,b)", title="Marcum Q-function")

```
## Parabolic cylinder functions

The parabolic cylinder functions `U(a, x)` and `V(a, x)` are fundamental solutions to the second-order linear differential equation
```math
    \frac{\text{d}^2 w}{\text{d}x^2} + \left( ax^2+bx+c \right)w = 0.
```

### `U(a, x)`
The function `U(a, x)` is defined via the integral representation
```math
    U(a,x) = \frac{\text{e}^{-\frac{1}{4}x^2}}{\Gamma\left(\frac{1}{2}+a\right)} \int_0^\infty t^{a-\frac{1}{2}} \text{e}^{-\frac{1}{2}t^2-xt} \, \text{d}t,
```
This representation is valid for all real values of `a` and `x`. 
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
The function `V(a, x)` is a second, linearly independent solution to the same differential equation satisfied by `U(a, x)`. The implementation of `V(a, x)` in `FewSpecialFunctions.jl` uses a combination of convergent series and asymptotic expansions adapted from standard references.

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
    D_n(x)= \frac{n}{x^n} \int_0^x \frac{t^n}{\text{e}^t-1} \, \text{d}t
```
And
```@example
using Plots, FewSpecialFunctions, LaTeXStrings # hide
ENV["GKSwstype"] = "100" # hide

plot_font = "Computer Modern" # hide
default(fontfamily=plot_font,linewidth=2.5, framestyle=:box, label=nothing, grid=true,palette=:tab10) # hide
x = range(0,stop=25,length=1000)
plot(x,Debye_function.(1,x),label=L"D_1(x)")
plot!(x,Debye_function.(2,x),label=L"D_2(x)")
plot!(x,Debye_function.(3,x), label=L"D_3(x)")
title!("Debye Functions")
xlabel!(L"x")
savefig("debye.svg"); nothing # hide
```
![](./debye.svg)




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
using Plots, FewSpecialFunctions, LaTeXStrings # hide
ENV["GKSwstype"] = "100" # hide

plot_font = "Computer Modern" # hide
default(fontfamily=plot_font,linewidth=2.5, framestyle=:box, label=nothing, grid=true,palette=:tab10) # hide
x = range(-5,stop=5,length=1000)
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
The Fresnel functions are both implemented using the trigonometric functions and the error function.

```math
    S(z) = \sqrt{\frac{\pi}{2}}\frac{1+i}{4} \bigg( \text{erf}\big(\frac{1+i}{\sqrt{2}}z \big) - i \text{erf}\big(\frac{1-i}{\sqrt{2}}z \big)\bigg)
```
And 
```math
    C(z) = \sqrt{\frac{\pi}{2}}\frac{1-i}{4} \bigg( \text{erf}\big(\frac{1+i}{\sqrt{2}}z \big) + i \text{erf}\big(\frac{1-i}{\sqrt{2}}z \big)\bigg)
```
The two implementations are shown in the examples below
```@example
using Plots, FewSpecialFunctions, LaTeXStrings # hide
ENV["GKSwstype"] = "100" # hide

plot_font = "Computer Modern" # hide
default(fontfamily=plot_font,linewidth=2.5, framestyle=:box, label=nothing, grid=true,palette=:tab10) # hide
x = range(-25,stop=25,length=5000)
plot(x,Fresnel_C_integral.(x),label=L"C(x)")
plot!(x,Fresnel_C_erf.(x), ls=:dash, lw=1.5, label=L"\tilde{C}(x)")
title!("Fresnel Integral")
xlabel!(L"x")
```
and

```@example
using Plots, FewSpecialFunctions, LaTeXStrings # hide
ENV["GKSwstype"] = "100" # hide

plot_font = "Computer Modern" # hide
default(fontfamily=plot_font,linewidth=2.5, framestyle=:box, label=nothing, grid=true,palette=:tab10) # hide
x = range(-25,stop=25,length=5000)
plot(x,Fresnel_S_integral.(x),label=L"S(x)")
plot!(x,Fresnel_S_erf.(x), ls=:dash, lw=1.5, label=L"\tilde{S}(x)")
title!("Fresnel Integral")
xlabel!(L"x")
```

## Hypergeometric functions

The confluent hypergeometric functions are solutions of Kummer’s equation
```math
    z \frac{d^2 F}{dz^2} +(b-z)\frac{dF}{dz}-a F = 0.
```
Kummer's equations has two linearly independent solutions given by
```math
    {}_1 F_{1}(a;b;z) 
```
and 
```math
    z^{1-b}{}_1F_1(a-b+1;2-b;z).
```
The numerical implementation in based on the following series expansion 
```math
    {}_1 F_1 (a,b,z) = 1+\frac{a}{b}z+\frac{a(a+1)}{b(b+1)2!}z^2 + \frac{a(a+1)(a+2)}{b(b+1)b(b+2)3!}z^3 + \cdots
```
Using this expansion the terms can be computed as
```math
    t(n) = \frac{(a+n-1)z}{(b+n-1)n}t(n-1),
```
and a cumulative product is produced. The errors are quite small; approximately ``10^{-12}-10^{-13}``. The figure below shows four examples of ``{}_1F_1``.
```@example
using Plots, FewSpecialFunctions, LaTeXStrings # hide
ENV["GKSwstype"] = "100" # hide

plot_font = "Computer Modern" # hide
default(fontfamily=plot_font,linewidth=2.5, framestyle=:box, label=nothing, grid=true,palette=:tab10) # hide
x = range(0,stop=3.5,length=100)
plot(x,confluent_hypergeometric_1F1.(2,0.5,x),label=L"{}_1F_1(2,0.5,x)")
plot!(x,confluent_hypergeometric_1F1.(3,0.7,x),label=L"{}_1F_1(3,0.7,x)")
plot!(x,confluent_hypergeometric_1F1.(2,1.1,x),label=L"{}_1F_1(2,1.1,x)")
plot!(x,confluent_hypergeometric_1F1.(3,1.3,x),label=L"{}_1F_1(3,1.3,x)")
title!("Confluent hypergeometric function")
xlabel!(L"x")
```

## Fermi-Dirac integrals

In solid state physics the Fermi-Dirac integral is given by

```math
    F_j(x) = \int_0^\infty \frac{t^j}{\exp(t-x)+1} \, dt.
```
Approximations to this and the normalized case for ``j=-1/2``, ``j=1/2``, ``j=3/2`` and ``j=5/2`` are implemented to varying accuacy. Most are of the order of ``10^{-12}``.

Here is an example

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
title!("Fermi-Dirac Integral")
xlabel!(L"x")
```
