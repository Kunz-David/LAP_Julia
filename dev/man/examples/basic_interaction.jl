
# # [Basic Interaction](@id basic_interaction)
#
#md # [![](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__/generated/example_placeholder.ipynb)
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/generated/example_placeholder.ipynb)

# In this example page I will show how to work with basic function that are meant for working with this registration algorithms.

# ## First how to generate some images that we can register
# We can do this with the [`gen_init`](@ref) for example.

## using the package
using LAP_julia;

## default arguments:
img, imgw, flow = gen_init();

## get a chess board image with a uniform flow:
img_chess, imgw_chess, flow_uniform = gen_init(:chess, :uniform, flow_args=[1 + 1im]);

# See also: [`gen_chess`](@ref), [`gen_quad_flow`](@ref), [`gen_quad_flow`](@ref), [`gen_lena`](@ref), [`gen_chess`](@ref), [`gen_anhir`](@ref)

# ## Next we would like to display these generated images and flows
# For that we will use the [`imgshow`](@ref), [`showflow`](@ref) and [`imgoverlay`](@ref) functions.

# See a single image:
imgshow(img)
#-----------
# See the default random flow:
showflow(flow)
#-----------
# See the uniform flow:
showflow(flow_uniform)
#-----------
# See the differences between 2 images:
imgoverlay(img, imgw)

# See also: [`imgshowflow`](@ref), [`warp_imgshowflow`](@ref), [`gen_quad_flow`](@ref), [`gen_lena`](@ref), [`gen_chess`](@ref), [`gen_anhir`](@ref),
