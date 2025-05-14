# FewSpecialFunctions

[![Documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://martinmikkelsen.github.io/FewSpecialFunctions.jl/dev/)
[![CI](https://github.com/MartinMikkelsen/FewSpecialFunctions.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/MartinMikkelsen/FewSpecialFunctions.jl/actions/workflows/ci.yml)
[![Coverage Status](https://coveralls.io/repos/github/MartinMikkelsen/FewSpecialFunctions.jl/badge.svg?branch=main)](https://coveralls.io/github/MartinMikkelsen/FewSpecialFunctions.jl?branch=main)
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)

A Julia package collecting a few special functions. Now includes over 13000 tests.

### Install 

Get the latest stable release with Julia's package manager:

```
julia ] add FewSpecialFunctions
```
Or use 
```
julia ] add https://github.com/MartinMikkelsen/FewSpecialFunctions.jl
```

### Examples
```julia
julia> using FewSpecialFunctions

julia> FermiDiracIntegral(3 / 2, 1.0)
2.6616826247307124
```

### Some other examples
![CombinedPlot](combinedplot.png)

### Cite

```yaml
cff-version: 1.2.0
title: FewSpecialFunctions.jl
message: >-
  You can cite if you want.
type: software
authors:
  - given-names: Martin
    family-names: Mikkelsen
    email: martin.mikkelsen@di.ku.dk
repository-code: 'https://github.com/MartinMikkelsen/FewSpecialFunctions.jl'
license: MIT
version: 0.1.3
date-released: '2023-01-26'
```
