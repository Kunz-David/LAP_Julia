```@meta
EditURL = "<unknown>/docs/src/man/examples/registration_functions.jl"
```

# Registration Functions Showcase

[![](https://mybinder.org/badge_logo.svg)](<unknown>/generated/registration_functions.ipynb)
[![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](<unknown>/generated/registration_functions.ipynb)

In this example page I will show how to work with the registration algorithms [`lap`](@ref), [`pflap`](@ref), [`sparse_lap`](@ref) and [`sparse_pflap`](@ref).

I will generate some images to show the workings and function calls.

```@example registration_functions
# using the package
using LAP_julia;

# default arguments:
img, imgw, flow = gen_init();
nothing #hide
```

See the [`Basic Interaction`](@ref basic_interaction) section or the [`Public Documentation`](@ref public_api) for more advanced input generation

## Full LAP Methods

### [`lap`](@ref) function

```@example registration_functions
# set a filter half size, which is larger than the default max displacement: 10 used in gen_init
fhs = 15;
flow_est, source_reg = lap(img, imgw, fhs);
nothing #hide
```

This gives us a estimation of the flow and the registered source image

```@eval

```@example registration_functions
pwd()
```

```

```@eval
using PyPlot
function optable(img1, img2, basename1, basename2, descr1, descr2)
    ## save fig 1 as next in line
    imgshow(img1); i = 1;
    fname1 = joinpath("..", "assets", string(basename1, i, ".png"))
    while isfile(fname1)
        i = i + 1
        fname1 = joinpath("..", "assets", string(basename1, i, ".png"))
    end
    savefig(fname1)
    ## save fig 2 as next in line
    imgshow(img2); j = 1;
    fname2 = joinpath("..", "assets", string(basename2, j, ".png"))
    while isfile(fname2)
        i = i + 1
        fname2 = joinpath("..", "assets", string(basename2, j, ".png"))
    end
    savefig(fname2)
    tbl = string(
        "$descr1 | $descr2\n",
        "------|--------\n",
        "![input]($fname1) | ![output]($fname2)\n"
    )
    Markdown.parse(tbl)
end

optable(img, imgw, "img", "imgw", "target", "source")
```

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

