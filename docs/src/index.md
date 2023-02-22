# FewSpecialFunctions.jl

*Some special functions.*

## Overview

[FewSpecialFunctions.jl](https://github.com/MartinMikkelsen/FewSpecialFunctions.jl) provides implementations of a few special functions. So far this includes
1. [Clausen functions](@ref)
2. [Debye functions](@ref)
3. [Regular Coulomb wave functions](@ref)
4. [Struve functions](@ref)
5. [Fresnel functions](@ref)

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
default()

x = range(-25,25,5000)

plot(Fresnel_C_erf.(x),Fresnel_S_erf.(x))
xlabel!(L"C(x)")
ylabel!(L"S(x)")
title!("Euler Spiral")
plot(Fresnel_C_erf.(x),Fresnel_S_erf.(x))
```
