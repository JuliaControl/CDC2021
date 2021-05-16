This repository contains the code that was submitted in our CDC2021 paper. The code should be run on the "dev" branch of ControlSystems.jl.

The following snippet will install the specific version of the ControlSystems.jl package and all packages used in the examples.
```julia
using Pkg
Pkg.add(PackageSpec(name="ControlSystems", rev="6853c47cdc4c87f381e1ef906a3449f3ad51e060"))
Pkg.add(["Plots", "GenericLinearAlgebra", "OrdinaryDiffEq", "DelayDiffEq", "Polynomials", "Symbolics", "CUDA", "BenchmarkTools", "SymPy", "MonteCarloMeasurements", "Optim"]),
```
Special care might be needed to get SymPy to work if you do not have SymPy installed in your python environment. See https://github.com/JuliaPy/SymPy.jl and related instructions.
