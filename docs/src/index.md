# FewSpecialFunctions.jl

*Some special functions.*

## Installation

Get the latest stable release with Julia's package manager:

```
julia ] add FewSpecialFunctions
```

## Features

- Numerical implementations of special mathematical functions
- **Generic type support**: all functions work with `Float32`, `Float64`, `BigFloat`, and other `AbstractFloat` types
- Automatic type promotion for mixed inputs

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

xlabel!(L"C(x)")
ylabel!(L"S(x)")
title!("Euler Spiral")
plot(FresnelC.(x),FresnelS.(x))
```
