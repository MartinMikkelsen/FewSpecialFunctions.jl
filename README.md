# FewSpecialFunctions

[![Documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://martinmikkelsen.github.io/FewSpecialFunctions.jl/dev/)
[![CI](https://github.com/MartinMikkelsen/FewSpecialFunctions.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/MartinMikkelsen/FewSpecialFunctions.jl/actions/workflows/ci.yml)
[![license-img]]: http://img.shields.io/badge/license-MIT-brightgreen.svg?style=flat


A Julia package collecting a few special functions.

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
    email: MartinMikkelsen@proton.me
repository-code: 'https://github.com/MartinMikkelsen/FewSpecialFunctions.jl'
license: MIT
version: 0.1.3
date-released: '2023-01-26'
```
