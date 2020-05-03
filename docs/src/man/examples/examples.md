# Examples

```@contents
Pages = ["single_lap.md", "polyfilter_lap.md"]
Depth = 5
```

Here I will show how the basic methods work on different kinds of displacements.
I use the word displacement and flow interchangeably. To see how to display the results check out the [Visualisation Functions](@ref) section.


## Generate test flows

### Unifrom
```@example uniform
using LAP_julia # hide
# Set some parameters:
ex_size = (200, 200)
ex_max_magnitude = 20
ex_tile_size = 999 # This makes the flow uniform

# Generate a random uniform flow
flow = gen_rand_flow(ex_size, ex_max_magnitude, ex_tile_size);

# Lets check what the flow looks like:
showflow(flow)
```

### Non-uniform Smooth
```@example smooth
using LAP_julia # hide
# Set some parameters:
ex_size = (200, 200)
ex_max_magnitude = 20
ex_tile_size = 60

# Generate a random uniform flow
flow = gen_rand_flow(ex_size, ex_max_magnitude, ex_tile_size);

# Lets check what the flow looks like:
showflow(flow)
```

_Check out the used functions: [`gen_rand_flow`](@ref), [`showflow`](@ref)_

## Generate test images

**Chessboard**
```@example uniform
chess = gen_chess()
```

```@example uniform
imgshow(chess)
```

**Lena image**
```@example uniform
using TestImages
lena = testimage("lena_gray")
```

```@setup smooth
using TestImages
chess = gen_chess()
lena = testimage("lena_gray")
```
