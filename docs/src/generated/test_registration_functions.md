```@meta
EditURL = "<unknown>/docs/src/man/examples/test_registration_functions.jl"
```

# Test Registration Functions

[![](https://mybinder.org/badge_logo.svg)](<unknown>/generated/test_registration_functions.ipynb)
[![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](<unknown>/generated/test_registration_functions.ipynb)

In this example page I will show how the registration algorithms ([`lap`](@ref), [`pflap`](@ref), [`sparse_lap`](@ref) and [`sparse_pflap`](@ref))
perform. I will will show their speed and the quality of their outputs.

First I some testing images:

```@example test_registration_functions
# using the package and TimerOutputs
using LAP_julia, TimerOutputs

# default arguments:
img, imgw, flow = gen_init();
nothing #hide
```

See the [`Basic Interaction`](@ref basic_interaction) section or the [`Public Documentation`](@ref public_api) for more advanced input generation

__These are the differences between the `target` (`img`) and `source` (`imgw`) images.__

```@example test_registration_functions
imgoverlay(img, imgw, figtitle="Target vs Source")
```

## [`test_registration_alg`](@ref)

This function times the registration algorithm on the inside and then test the output flow - `flow_est` and aligned source - `source_reg`.

```@example test_registration_functions
# choose a method
method = sparse_pflap
# start a timer
timer = TimerOutput("ALG: " * string(method))
# set the keyword arguments of the method
method_kwargs = Dict(:display => false, :timer => timer, :match_source_histogram => false)
# set the arguments of the method if any
method_args = []
# start the test
flow_est, source_reg, timer, results = test_registration_alg(method,
                                                             img,
                                                             imgw,
                                                             flow,
                                                             timer=timer,
                                                             method_args=method_args,
                                                             only_flow_compare=false, # this adds source reg tests.
                                                             method_kwargs=method_kwargs);
nothing #hide
```

#### Check the generated flow and the original flow:
__Ground Truth:__

```@example test_registration_functions
showflow(flow.*(-1), figtitle="Ground Truth")
```

__Estimated Flow:__

```@example test_registration_functions
showflow(flow_est, figtitle="Estimated Flow")
```

### Compare it with the [`pflap`](@ref) function

```@example test_registration_functions
method = pflap
# same inputs, new timer
timer = TimerOutput("ALG: " * string(method))
flow_est, source_reg, timer, results = test_registration_alg(method,
                                                             img,
                                                             imgw,
                                                             flow,
                                                             timer=timer,
                                                             method_args=method_args,
                                                             only_flow_compare=false, # this adds source reg tests.
                                                             method_kwargs=method_kwargs);
nothing #hide
```

#### Check the generated flow and the original flow:
__Ground Truth:__

```@example test_registration_functions
showflow(flow.*(-1), figtitle="Ground Truth")
```

__Estimated Flow:__

```@example test_registration_functions
showflow(flow_est, figtitle="Estimated Flow")
```

The [`time_reg_alg`](@ref) is used to time the registration algorithm, when the original flow isn't known.

## using the `display` keyword
The registration algorithms have a display keyword argument, that can be set to true to output figures
showing the `flow_est` at each iteration and print some extra debug info.

```@example test_registration_functions
method = sparse_pflap
timer = TimerOutput("ALG: " * string(method))
method_kwargs = Dict(:display => true, :timer => timer, :match_source_histogram => false)
method_args = []
flow_est, source_reg, timer, results, (figs,) = test_registration_alg(method,
                                                                      img,
                                                                      imgw,
                                                                      flow,
                                                                      timer=timer,
                                                                      method_args=method_args,
                                                                      only_flow_compare=false, # this adds source reg tests.
                                                                      method_kwargs=method_kwargs);
nothing #hide
```

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

