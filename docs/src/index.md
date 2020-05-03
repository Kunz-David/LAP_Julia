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

```@example
a = 1
b = 2
a + b
```


```@example
using LAP_julia
flow = gen_rand_flow()
showflow(flow)
```

```julia
# some julia code
d = 12
```

```@contents
```

## Usage

## Index

```@index
```
