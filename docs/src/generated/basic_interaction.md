```@meta
EditURL = "<unknown>/docs/src/man/examples/basic_interaction.jl"
```

# Basic interaction

[![](https://mybinder.org/badge_logo.svg)](<unknown>/generated/basic_interaction.ipynb)
[![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](<unknown>/generated/basic_interaction.ipynb)

In this example page I will show how to work with basic function that are meant for working with this registration algorithms.

## First how to generate some images that we can register
We can do this with the [`gen_init`](@ref) for example.

```@example basic_interaction
# using the package
using LAP_julia;

# default arguments:
img, imgw, flow = gen_init();

# get a chess board image with a uniform flow:
img_chess, imgw_chess, flow_uniform = gen_init(:chess, :uniform, flow_args=[1 + 1im]);
nothing #hide
```

See also: [`gen_chess`](@ref), [`gen_quad_flow`](@ref), [`gen_quad_flow`](@ref), [`gen_lena`](@ref), [`gen_chess`](@ref), [`gen_anhir`](@ref)

## Next we would like to display these generated images and flows
For that we will use the [`imgshow`](@ref), [`showflow`](@ref) and [`imgoverlay`](@ref) functions.

# to see a single image

```@example basic_interaction
imgshow(img)
```

# to see the default random flow

```@example basic_interaction
showflow(flow)
```

# to see the uniform flow

```@example basic_interaction
showflow(flow_uniform)
```

# to see the differences between 2 images

```@example basic_interaction
imgoverlay(img, imgw)
```

See also: [`imgshowflow`](@ref), [`warp_imgshowflow`](@ref), [`gen_quad_flow`](@ref), [`gen_lena`](@ref), [`gen_chess`](@ref), [`gen_anhir`](@ref),

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

