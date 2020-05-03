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

```@contents
Pages = [
    "man/guide.md",
    "man/examples.md",
    "man/syntax.md",
    "man/doctests.md",
    "man/hosting.md",
    "man/latex.md",
    "man/contributing.md",
]
Depth = 1
```


## Library Outline

```@contents
Pages = ["lib/public.md", "lib/private.md"]
```

<!-- ```@example
using LAP_julia
flow = gen_rand_flow()
showflow(flow)
``` -->

```@contents
```
