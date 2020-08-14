# LAP_julia.jl

__Image registration in Julia__

## Installation

To install paste this into a Julia terminal:
```julia
using Pkg; Pkg.add(PackageSpec(url="https://github.com/Kunz-David/LAP_julia"))
```

!!! note "Linux"

    The plotting functions use Julia's PyPlot module and Matplotlib has to be installed in your default Python. You can either install Matplotlib in your Python or let Julia use the Python it installed and manages. For that set ENV[PYTHON] to the Python installed by Julia. So something like this, will do the trick:

    ```julia
    ENV["PYTHON"] = "..ENTER USER.../.julia/conda/3/bin"
    ```

## Manual Outline

To see this section, visit the [online docs](https://kunz-david.github.io/LAP_Julia.jl/dev/), where there are examples of how to run the methods.


## Library Outline

```@contents
Pages = ["lib/public.md", "lib/private.md"]
```
