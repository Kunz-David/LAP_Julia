# Examples

```@contents
pages = ["examples.md"]
Depth = 3
```

Here I will show how the basic methods work on different kinds of displacements.
I use the word displacement and flow interchangeably.


## Generate test flows

### Unifrom
```@example uniform
using LAP_julia # hide
using PyPlot
# Set some parameters:
ex_size = (256, 256)
ex_max_magnitude = 20
ex_tile_size = 999 # This makes the flow uniform

# Generate a random uniform flow
uni_flow = gen_tiled_flow(ex_size, ex_max_magnitude, ex_tile_size);

# Lets check what the flow looks like:
showflow(uni_flow, figtitle="Uniform warp")
```

*Check out the used functions: [`gen_tiled_flow`](@ref), [`showflow`](@ref)*

### Non-uniform Smooth
```@example smooth
using LAP_julia # hide
using PyPlot
# Set some parameters:
ex_size = (256, 256)
ex_max_magnitude = 20
ex_tile_size = 60

# Generate a random uniform flow
smooth_flow = gen_tiled_flow(ex_size, ex_max_magnitude, ex_tile_size);

# Lets check what the flow looks like:
showflow(smooth_flow, figtitle="Smooth warp")
```

*Check out the used functions: [`gen_tiled_flow`](@ref), [`showflow`](@ref)*

## Generate test images

### Chessboard
**Normal**
```@example uniform
chess = gen_chess(64, 4)
imgshow(chess, figtitle="Chessboard")
```
*Check out the used functions: [`gen_chess`](@ref)*

```@setup smooth
chess = gen_chess(64, 4)
```

**Warped with `uni_flow`**
```@example uniform
chess_uni_warped = warp_img(chess, -real(uni_flow), -imag(uni_flow))
imgshowflow(chess_uni_warped, uni_flow, figtitle="Chessboard uniform warp")
```

**Warped with `smooth_flow`**
```@example smooth
chess_smooth_warped = warp_img(chess, -real(smooth_flow), -imag(smooth_flow))
imgshowflow(chess_smooth_warped, smooth_flow, figtitle="Chessboard smooth warp")
```
*Check out the used functions: [`warp_img`](@ref), [`imgshowflow`](@ref)*

### Lena
**Normal**
```@example uniform
using TestImages
lena = testimage("lena_gray")
# The Lena image from TestImages has elements of type ColorTypes.Gray{FixedPointNumbers.Normed{UInt8,8}}
# so to use our algorithms later we need to convert them to floats.
lena = Float32.(lena)

imgshow(lena, figtitle="Lena")
```

```@setup smooth
using TestImages
lena = testimage("lena_gray")

# The Lena image from TestImages has elements of type ColorTypes.Gray{FixedPointNumbers.Normed{UInt8,8}} so we convert it to Float32.
lena = Float32.(lena)
```

**Warped with `uni_flow`**
```@example uniform
lena_uni_warped = warp_img(lena, -real(uni_flow), -imag(uni_flow))
imgshowflow(lena_uni_warped, uni_flow, figtitle="Lena uniform warp")
```

**Warped with `smooth_flow`**
```@example smooth
lena_smooth_warped = warp_img(lena, -real(smooth_flow), -imag(smooth_flow))
imgshowflow(lena_smooth_warped, smooth_flow, figtitle="Lena smooth warp")
```
*Check out the used functions: [`warp_img`](@ref), [`imgshowflow`](@ref)*

-------







# `single_lap` function

Here are examples of using the [`single_lap`](@ref) function to find the flow between the original and the warped image.
I will show the results of the algorithm for an image warped by a uniform flow and a smoothly varying flow. (The code running is the same.)

## Uniform Flow

This how the algorithm performs on a uniform flow.
First we have to choose a `filter_half_size` which has to be the same or higher that the `ex_max_magnitude` chosen.
Then we choose the `window_size` parameter of the algorithm, usually a list like this: `[2*filter_half_size + 1, 2*filter_half_size + 1]`.

```@setup uniform
# save original
showflow(uni_flow, figtitle="Original Uniform Flow")
savefig("orig_uni_flow.png")
```

### Chessboard
Run:
```@example uniform
filter_half_size = 20
window_size = [2*filter_half_size + 1, 2*filter_half_size + 1]

estim_flow = single_lap(chess, chess_uni_warped, filter_half_size, window_size)
nothing # hide
```
*Check out [`single_lap`](@ref).*

Check the results:
```@example uniform
# save estimation
showflow(estim_flow, figtitle="Single Estimated Uniform Flow")
savefig("single_chess_estim_uni_flow.png")
nothing # hide
```

| Original        | Estimated     |
|:---------------:|:-------------:|
| ![](orig_uni_flow.png) | ![](single_chess_estim_uni_flow.png) |

### Lena

Run:
```@example uniform
filter_half_size = 20
window_size = [2*filter_half_size + 1, 2*filter_half_size + 1]

estim_flow = single_lap(lena, lena_uni_warped, filter_half_size, window_size)
nothing # hide
```
*Check out [`single_lap`](@ref).*

Check the results:
```@example uniform
# save estimation
showflow(estim_flow, figtitle="Single Estimated Uniform Flow")
savefig("single_lena_estim_uni_flow.png")
nothing # hide
```

| Original        | Estimated     |
|:---------------:|:-------------:|
| ![](orig_uni_flow.png) | ![](single_lena_estim_uni_flow.png) |





## Non-uniform Smooth Flow

This how the algorithm performs on a non-uniform smoothly varying flow.

```@setup smooth
# save original
showflow(smooth_flow, figtitle="Original smooth flow")
savefig("orig_smooth_flow.png")
```

### Chessboard
Run:
```@example smooth
filter_half_size = 20
window_size = [2*filter_half_size + 1, 2*filter_half_size + 1]

estim_flow = single_lap(chess, chess_smooth_warped, filter_half_size, window_size)
nothing # hide
```

Check the results:
```@example smooth
# save estimation
showflow(estim_flow, figtitle="Single Estimated Smooth Flow")
savefig("single_chess_estim_smooth_flow.png")
nothing # hide
```

| Original        | Estimated     |
|:---------------:|:-------------:|
| ![](orig_smooth_flow.png) | ![](single_chess_estim_smooth_flow.png) |

### Lena

Run:
```@example smooth
filter_half_size = 20
window_size = [2*filter_half_size + 1, 2*filter_half_size + 1]

estim_flow = single_lap(lena, lena_smooth_warped, filter_half_size, window_size)
nothing # hide
```

Check the results:
```@example smooth
# save estimation
showflow(estim_flow, figtitle="Single Estimated Smooth Flow")
savefig("single_lena_estim_smooth_flow.png")
nothing # hide
```

| Original        | Estimated     |
|:---------------:|:-------------:|
| ![](orig_smooth_flow.png) | ![](single_lena_estim_smooth_flow.png) |

-----------



# `polyfilter_lap` function

Here are examples of using the [`polyfilter_lap`](@ref) function to find the flow between the original and the warped image.
I will show the results of the algorithm for an image warped by a uniform flow and a smoothly varying flow. (The code running is the same.)

## Uniform Flow

This how the algorithm performs on a uniform flow.

### Chessboard
Run:
```@example uniform
filter_half_size = 20
window_size = [2*filter_half_size + 1, 2*filter_half_size + 1]

estim_flow, source_reg = polyfilter_lap(chess, chess_uni_warped, display=false)
nothing # hide
```
*Check out [`polyfilter_lap`](@ref).*

Check the results:
```@example uniform
# save estimation
showflow(estim_flow, figtitle="Polyfilter Estimated Uniform Flow")
savefig("polyfilter_chess_estim_uni_flow.png")
nothing # hide
```

| Original        | Estimated     |
|:---------------:|:-------------:|
| ![](orig_uni_flow.png) | ![](polyfilter_chess_estim_uni_flow.png) |

### Lena

Run:
```@example uniform
filter_half_size = 20
window_size = [2*filter_half_size + 1, 2*filter_half_size + 1]

estim_flow, source_reg = polyfilter_lap(lena, lena_uni_warped, display=false)
nothing # hide
```
*Check out [`polyfilter_lap`](@ref).*

Check the results:
```@example uniform
# save estimation
showflow(estim_flow, figtitle="Polyfilter Estimated Uniform Flow")
savefig("polyfilter_lena_estim_uni_flow.png")
nothing # hide
```

| Original        | Estimated     |
|:---------------:|:-------------:|
| ![](orig_uni_flow.png) | ![](polyfilter_lena_estim_uni_flow.png) |





## Non-uniform Smooth

This how the algorithm performs on a non-uniform smooth flow.

### Chessboard
Run:
```@example smooth
estim_flow, source_reg = polyfilter_lap(chess, chess_smooth_warped, display=false)
nothing # hide
```

Check the results:
```@example smooth
# save estimation
showflow(estim_flow, figtitle="Polyfilter Estimated Smooth Flow")
savefig("polyfilter_chess_estim_smooth_flow.png")
nothing # hide
```

| Original        | Estimated     |
|:---------------:|:-------------:|
| ![](orig_smooth_flow.png) | ![](polyfilter_chess_estim_smooth_flow.png) |

### Lena

Run:
```@example smooth
filter_half_size = 20
window_size = [2*filter_half_size + 1, 2*filter_half_size + 1]

estim_flow, source_reg = polyfilter_lap(lena, lena_smooth_warped, display=false)
nothing # hide
```

Check the results:
```@example smooth
# save estimation
showflow(estim_flow, figtitle="Polyfilter Estimated Smooth Flow")
savefig("polyfilter_lena_estim_smooth_flow.png")
nothing # hide
```

| Original        | Estimated     |
|:---------------:|:-------------:|
| ![](orig_smooth_flow.png) | ![](polyfilter_lena_estim_smooth_flow.png) |
