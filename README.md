# FewSpecialFunctions

[![Documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://martinmikkelsen.github.io/FewSpecialFunctions.jl/dev/)


Implementation of a few special functions. So far this includes

1. Clausen function 
2. Coulomb wave functions
3. Debye function
4. Fresnel functions
5. Struve functions
6. Hypergeometric functions
7. Confluent Hypergeometric functions
8. Fermi-Dirac Integrals

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

![CombinedPlot](combinedplot.png)

### Contribute

If you find a problem with the functions or docs, please open an issue. If you want to implement a function please make a PR.

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