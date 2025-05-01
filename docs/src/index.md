# FewSpecialFunctions.jl

*Some special functions.*

## Installation

Get the latest stable release with Julia's package manager:

```
julia ] add FewSpecialFunctions
```

## Quick example

Here is how to generate an Euler spiral using [FewSpecialFunctions.jl](https://github.com/MartinMikkelsen/FewSpecialFunctions.jl). 

```@example EulerSpiral
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

x = range(-25,25,length=5000)

plot(Fresnel_C_erf.(x),Fresnel_S_erf.(x))
xlabel!(L"C(x)")
ylabel!(L"S(x)")
title!("Euler Spiral")
plot(Fresnel_C_erf.(x),Fresnel_S_erf.(x))
```
