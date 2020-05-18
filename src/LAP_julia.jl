
module LAP_julia

export
    # visualisation funcitons
    showflow, imgshowflow, imgshow, warp_imgshowflow, showsparseflow,
    # useful helpers
    mean,
    # data generation
    gen_rand_flow, gen_chess,
    # Main LAP funcitons
    single_lap, polyfilter_lap,
    polyfilter_lap_at_points, single_lap_at_points,
    # Types
    SimpleKeypoint, Image, Flow,
    # Point location
    find_edge_points,
    # Interpolation
    warp_img, interpolate_flow


"""
    Image{T} = Matrix{T} where T <: Real

Image is a `Matrix` with elements that are `Real`.
"""
const Image{T} = Matrix{T} where T <: Real

"""
    Flow{T} = Matrix{Complex{T}} where T <: Real

Flow is a `Matrix` with elements that are `Complex`.
"""
const Flow{T} = Matrix{Complex{T}} where T <: Real


include("helpers.jl")
include("inpaint.jl")
include("lap.jl")
include("interpolation.jl")
include("visualise.jl")
include("data_gen.jl")
include("gradient_points.jl")

# Revise.includet("lap.jl")

using .lap
using .visualise
using .inpaint: inpaint_nans!
using .gradient_points
# using .visualise: showflow, imgshowflow, imgshow
using .data_gen
using .interpolation


loaded() = println("LAP_julia succesfully loaded!")

# some macros:
# using BenchmarkTools
# macro debugtime(body, note)
#     quote
#         println($note)
#         @btime $body
#     end
# end
# # @debugtime(1+1, "1 + 1")
#
# macro debugassert(body, note)
#     quote
#         println($note)
#         @assert $body
#     end
# end
# # @debugassert (1+1 == 1) "math works"

loaded()


end # module
