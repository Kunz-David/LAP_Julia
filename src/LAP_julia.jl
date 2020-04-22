
module LAP_julia

export showflow, gen_rand_flow, img_showflow, imgshow, mean
export single_lap, polyfilter_lap, SimpleKeypoint, find_keypoints_from_gradients

struct SimpleKeypoint
    pos::Tuple{T, T} where T <: Int
    value::T where T <: Number
end

Base.show(io::IO, m::SimpleKeypoint) = print(io, m.pos, " ", m.value, )


include("helpers.jl")
include("lap.jl")
include("inpaint.jl")
include("visualise.jl")
include("data_gen.jl")
include("interpolation.jl")
include("gradient_points.jl")

# Revise.includet("lap.jl")

using .lap: single_lap, polyfilter_lap
# using .inpaint
# using .helpers
using .gradient_points
using .helpers: mean
using .visualise: showflow, img_showflow, imgshow
using .data_gen: gen_rand_flow
# using .interpolation


greet() = println("Hello             World!")

# some macros:

using BenchmarkTools
macro debugtime(body, note)
    quote
        println($note)
        @btime $body
    end
end
# @debugtime(1+1, "1 + 1")

macro debugassert(body, note)
    quote
        println($note)
        @assert $body
    end
end
# @debugassert (1+1 == 1) "math works"




# showall(x) = show(stdout, "text/plain", x)
# export showall

greet()


end # module
