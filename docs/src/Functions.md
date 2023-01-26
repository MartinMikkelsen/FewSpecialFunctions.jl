# Functions

## Clausen functions
The Clausen function is given by

``Cl_2(\phi)=-\int_0^\phi \log|2\sin(x/2)| dx``

# Examples 

```@example overview
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
plot(x,Clausen.(x), label=L"Cl_2(ϕ)")
xlabel!(L"ϕ")
title!("Clausen function")
```

## Debye functions

```@example overview
using Plots, FewSpecialFunctions

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
```

## Regular Coulomb wave functions
This implementation need the gamma function from [SpecialFunctions.jl](https://github.com/JuliaMath/SpecialFunctions.jl)
```@example overview

using Plots, FewSpecialFunctions

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
plot!(x,regular_coulomb.(0,-0.3,x), label=L"F_0(0.3,ρ)")
xlabel!(L"ρ")
title!("Regular Coulomb Wave Functions")
```
And 

```@example overview

using Plots, FewSpecialFunctions

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

```@example overview

using Plots, FewSpecialFunctions

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

plot(x,Struve.(0,x),label=L"H_0(x)")
plot!(x,Struve.(1,x),label=L"H_1(x)")
plot!(x,Struve.(2,x),label=L"H_3(x)")
xlabel!(L"x")
title!("Struve Functions")
```
## Fresnel functions

```@example overview

using Plots, FewSpecialFunctions

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
plot!(x,Fresnel_C_err.(x), ls=:dash, lw=1.5, label=L"\tilde{C}(x)")
title!("Fresnel Integral")
xlabel!(L"x")

```
and

```@example overview

using Plots, FewSpecialFunctions

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
plot!(x,Fresnel_S_err.(x), ls=:dash, lw=1.5, label=L"\tilde{S}(x)")
title!("Fresnel Integral")
xlabel!(L"x")

```
