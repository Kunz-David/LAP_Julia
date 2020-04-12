
module LAP_julia

include("helpers.jl")
include("lap.jl")
include("inpaint.jl")
include("visualise.jl")
include("data_gen.jl")
include("interpolation.jl")

# Revise.includet("lap.jl")

using .lap: single_lap, polyfilter_lap
using .inpaint
using .helpers
using .helpers: mean
using .visualise: showflow, img_showflow, imgshow
using .data_gen: gen_rand_flow
using .interpolation
# used public modules:
using ImageFiltering: imfilter!, kernelfactors, centered

export imshow, showflow, gen_rand_flow, img_showflow, imgshow, mean
export single_lap, polyfilter_lap

greet() = println("Hello             World!")

# showall(x) = show(stdout, "text/plain", x)
# export showall

greet()


end # module
