
# # Registration Functions
#
#md # [![](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__/generated/example_placeholder.ipynb)
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/generated/example_placeholder.ipynb)

# In this example page I will show how to work with the registration algorithms [`lap`](@ref), [`pflap`](@ref), [`sparse_lap`](@ref) and [`sparse_pflap`](@ref).

# I will generate some images to show the workings and function calls.

## using the package
using LAP_julia;

## default arguments:
img, imgw, flow = gen_init();

# See the [`Basic Interaction`](@ref basic_interaction) section or the [`Public Documentation`](@ref public_api) for more advanced input generation

# __These are the differences between the `target` (`img`) and `source` (`imgw`) images.__

imgoverlay(img, imgw, figtitle="Target vs Source")

# ## Full LAP Methods

# These functions are from the algorithms from the paper [Local All-Pass Geometric Deformations](https://www.researchgate.net/publication/320574179_Local_All-Pass_Geometric_Deformations).

# ### [`lap`](@ref) function

## set a filter half size, which is larger than the default max displacement: 10 used in gen_init
fhs = 15;
flow_est, source_reg = lap(img, imgw, fhs);

# This gives us a estimation of the flow and the registered source image

# __This is what the registered source image looks like:__
imgoverlay(img, source_reg, figtitle="lap: Target vs Registered Source")

# ### [`pflap`](@ref) function

flow_est, source_reg = pflap(img, imgw, display=false);

# This gives us a estimation of the flow and the registered source image

# __This is what the registered source image looks like:__
imgoverlay(img, source_reg, figtitle="pflap: Target vs Registered Source")

# ## Sparse LAP Methods

# These functions are insired by the [`lap`](@ref) and [`pflap`](@ref) functions.
# The difference being that the displacement vectors are calculated at points of high gradient and
# the result is fit into a global deformation.

# ### [`sparse_lap`](@ref) function

## set a filter half size, which is larger than the default max displacement: 10 used in gen_init
fhs = 15;
flow_est, source_reg = sparse_lap(img, imgw, fhs, display=false);

# This gives us a estimation of the flow and the registered source image

# __This is what the registered source image looks like:__
imgoverlay(img, source_reg, figtitle="sparse lap: Target vs Registered Source")

# ### [`sparse_pflap`](@ref) function

flow_est, source_reg = sparse_pflap(img, imgw, display=false);

# This gives us a estimation of the flow and the registered source image

# __This is what the registered source image looks like:__
imgoverlay(img, source_reg, figtitle="sparse pflap: Target vs Registered Source")
