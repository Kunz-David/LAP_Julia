
module LAP_julia

using TimerOutputs

export
    # visualisation funcitons
    showflow,
    imgshowflow,
    imgshow,
    warp_imgshowflow,
    showsparseflow,
    addpoints,
    imgoverlay,
    # useful helpers
    mean,
    classic_alg,
    create_sparse_flow_from_sparse,
    create_sparse_flow_from_full,
    # data generation
    gen_quad_flow,
    gen_tiled_flow,
    gen_uniform_flow,
    gen_chess,
    gen_init,
    gen_lena,
    gen_anhir,
    load_anhir_image_pair,
    # Main LAP funcitons
    lap,
    sparse_lap,
    sparse_lap_win_sum1,
    sparse_pflap,
    single_lap,
    pflap,
    sparse_pflap,
    single_lap_at_points,
    # Types
    Image,
    Flow,
    # Point location
    find_edge_points,
    # Interpolation
    warp_img,
    interpolate_flow,
    interpolate_flow_quad,
    #experimenting
    test_registration_alg,
    assess_flow_quality,
    assess_source_reg_quality


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
include("interpolation.jl")
include("visualise.jl")
include("data_gen.jl")
include("gradient_points.jl")
include("lap_algs/helpers.jl")
include("lap_algs/multi.jl")
include("lap_algs/single.jl")
include("experimenting/speedtest.jl")
include("experimenting/qualitytest.jl")
include("experimenting/visualise.jl")
include("birl.jl")


loaded() = println("LAP_julia succesfully loaded!")

loaded()


end # module
