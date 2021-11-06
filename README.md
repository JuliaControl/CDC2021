This repository contains the code that was submitted in our CDC2021 paper. The code should be run on the "dev" branch of ControlSystems.jl.

The following snippet will replicate the environment that was used to generate all figures if run from the repos base directory.
```
julia --project=code -E "import Pkg; Pkg.instantiate()"
```
Special care might be needed to get SymPy to work if you do not have SymPy installed in your python environment. See https://github.com/JuliaPy/SymPy.jl and related instructions.

To later use this environment to run one of the scripts you would do
```
julia --project=code code/somefile.jl
```