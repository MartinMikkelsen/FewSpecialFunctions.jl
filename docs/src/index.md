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

x = range(-25,25,5000)

plot(Fresnel_C_err.(x),Fresnel_S_err.(x))
xlabel!(L"C(x)")
ylabel!(L"S(x)")
title!("Euler Spiral")
```

## About

MIT License

Copyright (c) 2023 Martin Mikkelsen

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.