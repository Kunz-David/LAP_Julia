# Public Documentation

Documentation for `LAP_julia.jl`'s public interface.

See the Internals section of the manual for internal package docs covering all submodules.

## Contents

```@contents
Pages = ["public.md"]
```

## Index

```@index
Pages = ["public.md"]
```

## Public Interface
### Algorithm Functions
```@docs
single_lap
polyfilter_lap
```
### Visualisation Functions
```@docs
showflow
imgshowflow
imgshow
warp_imgshowflow
```
### Generate-Data Functions
```@docs
gen_rand_flow
gen_chess
```
### Types
```@docs
Image
Flow
```
### Interpolation
```@docs
warp_img
```

### Other
```@docs
LAP_julia.gradient_points.find_keypoints_from_gradients
```
