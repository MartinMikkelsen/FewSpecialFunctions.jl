# FewSpecialFunctions

[![Documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://martinmikkelsen.github.io/FewSpecialFunctions.jl/dev/)
[![CI](https://github.com/MartinMikkelsen/FewSpecialFunctions.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/MartinMikkelsen/FewSpecialFunctions.jl/actions/workflows/ci.yml)
[![codecov](https://codecov.io/gh/MartinMikkelsen/FewSpecialFunctions.jl/graph/badge.svg?token=M6UNBWGD18)](https://codecov.io/gh/MartinMikkelsen/FewSpecialFunctions.jl)
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)

A Julia package collecting a few special functions. Now includes over 13000 tests. Includes the following functions:

- [Clausen functions](https://en.wikipedia.org/wiki/Clausen_function)
- [Marcum-Q functions](https://en.wikipedia.org/wiki/Marcum_Q-function)
- [Parabolic cylinder functions](https://en.wikipedia.org/wiki/Parabolic_cylinder_function)
- [Coulomb wave functions](https://en.wikipedia.org/wiki/Coulomb_wave_function)
- [Debye functions](https://en.wikipedia.org/wiki/Debye_function)
- [Fermi-Dirac integrals](https://en.wikipedia.org/wiki/Incomplete_Fermi%E2%80%93Dirac_integral)
- [Fresnel integrals](https://en.wikipedia.org/wiki/Fresnel_integral)

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

Some other examples are shown in the [documentation](https://martinmikkelsen.github.io/FewSpecialFunctions.jl/dev/).

![CombinedPlot](combinedplot.png)

