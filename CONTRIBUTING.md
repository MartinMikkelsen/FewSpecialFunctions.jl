# Contributing to FewSpecialFunctions.jl

Bug reports, corrections, and new functions are welcome. This file describes
how to work on the package locally.

## Setup

The package supports Julia 1.10 and later. Clone the repository and install
its dependencies:

```
git clone https://github.com/MartinMikkelsen/FewSpecialFunctions.jl
cd FewSpecialFunctions.jl
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

## Running the tests

Run the full test suite, which includes the Aqua.jl quality checks and the
ForwardDiff extension tests:

```
julia --project=. -e 'using Pkg; Pkg.test()'
```

To run the tests of one function family, start `julia --project=.` in an
environment where the test dependencies (`Test`, `DelimitedFiles`,
`ForwardDiff`, `Aqua`) are available and include the file:

```julia
using FewSpecialFunctions, Test
include("test/test_Coulomb.jl")
```

New functions need tests against independent reference values (another
library, a closed form, or high-precision quadrature). State in the test or
the docstring where the reference values come from.

## Building the documentation

The manual is built with Documenter.jl from the `docs` environment:

```
julia --project=docs -e 'using Pkg; Pkg.instantiate()'
julia --project=docs docs/make.jl
```

The output is written to `docs/build`; open `docs/build/index.html` in a
browser. The build runs every `@example` block and every `jldoctest` block in
`docs/src` and in the docstrings, and fails if a doctest's output differs
from what is written. Plain `julia` code fences are not executed, so prefer
`jldoctest` for examples whose output matters. The deployment step at the end
of `docs/make.jl` does nothing outside continuous integration.

The error figures quoted in the manual and docstrings are printed by
`docs/accuracy_checks.jl`, which compares the package against closed forms,
series, and high-precision quadrature:

```
julia --project=docs docs/accuracy_checks.jl
```

If a change alters one of those figures, update the text to the new output.
When you quote a new error figure, add its measurement to that script and
state the test point or grid next to the figure.

Every exported function needs a docstring and an entry in
`docs/src/API.md`. A docstring should state the mathematical definition, the
accepted arguments and their domain, the errors thrown, and what is known
about the accuracy, distinguishing relative from absolute error.

## Formatting

The code is formatted with [Runic.jl](https://github.com/fredrikekre/Runic.jl).
Install it once into a shared environment and format the sources in place:

```
julia --project=@runic -e 'using Pkg; Pkg.add("Runic")'
julia --project=@runic -m Runic --inplace src ext test docs/make.jl
```

Use `--check --diff` instead of `--inplace` to see what would change. The
`julia -m` form requires Julia 1.12 or later; on older versions run
`julia --project=@runic -e 'using Runic; exit(Runic.main(ARGS))' -- --inplace src ext test docs/make.jl`.

## Reporting a bug

Open an issue at
<https://github.com/MartinMikkelsen/FewSpecialFunctions.jl/issues> with

- the smallest call that shows the problem, including the argument types
  (for example `FermiDiracIntegral(1, big"1.0")`),
- the value or error you obtained and the value you expected, with the source
  of the expected value,
- the output of `versioninfo()` and of `import Pkg; Pkg.status("FewSpecialFunctions")`.

For an accuracy problem, say whether the discrepancy is a relative or an
absolute error.

## Pull requests

Keep a pull request to one topic. Before opening it, run the tests, build the
documentation if you changed docstrings or `docs/`, and format the code.

## License

Contributions are distributed under the [MIT license](LICENSE) of the
package.
