module lap

using TimerOutputs

export polyfilter_lap, single_lap,
       polyfilter_lap_at_points, single_lap_at_points

include("lap_algs/helpers.jl")
include("lap_algs/single.jl")
include("lap_algs/multi.jl")

end # module
