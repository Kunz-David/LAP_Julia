
# # Test Registration Functions Showcase
#
#md # [![](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__/generated/example_placeholder.ipynb)
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/generated/example_placeholder.ipynb)

# In this example page I will show how the registration algorithms ([`lap`](@ref), [`pflap`](@ref), [`sparse_lap`](@ref) and [`sparse_pflap`](@ref))
# perform. I will will show their speed and the quality of their outputs.

# First I some testing images:

## using the package
using LAP_julia;

## default arguments:
img, imgw, flow = gen_init();

# See the [`Basic Interaction`](@ref basic_interaction) section or the [`Public Documentation`](@ref public_api) for more advanced input generation

# __These are the differences between the `target` (`img`) and `source` (`imgw`) images.__

imgoverlay(img, imgw, figtitle="Target vs Source")

# ## [`test_registration_alg`](@ref)

# This function times the registration algorithm on the inside and then test the output flow - `flow_est` and aligned source - `source_reg`.

## choose a method
method = sparse_pflap
## start a timer
timer = TimerOutput("ALG: " * string(method))
## set the keyword arguments of the method
method_kwargs = Dict(:display => false, :timer => timer, :match_source_histogram => false)
## set the arguments of the method if any
method_args = []
## start the test
flow_est, source_reg, timer, results = test_registration_alg(method,
                                                             img,
                                                             imgw,
                                                             flow,
                                                             timer=timer,
                                                             method_args=method_args,
                                                             only_flow_compare=false, # this adds source reg tests.
                                                             method_kwargs=method_kwargs)


# __Check the generated flow and the original flow:__
showflow(flow.*(-1))
showflow(flow_est)

# #### Compare it with the [`pflap`](@ref) function

method = pflap
## same inputs, new timer
timer = TimerOutput("ALG: " * string(method))
flow_est, source_reg, timer, results = test_registration_alg(method,
                                                             img,
                                                             imgw,
                                                             flow,
                                                             timer=timer,
                                                             method_args=method_args,
                                                             only_flow_compare=false, # this adds source reg tests.
                                                             method_kwargs=method_kwargs)


# __Check the generated flow and the original flow:__
showflow(flow.*(-1))
showflow(flow_est)

# The [`time_reg_alg`](@ref) is used to time the registration algorithm, when the original flow isn't known.


# ## using the `display` keyword
# The registration algorithms have a display keyword argument, that can be set to true to output figures
# showing the `flow_est` at each iteration and print some extra debug info.
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
                                                                      method_kwargs=method_kwargs)
