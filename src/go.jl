
using PyPlot
try
    using LAP_julia
catch e
    @warn e "fail 1st time"
    try
        using LAP_julia
    catch e
        @warn e "fail 2nd time"
        try
            using LAP_julia
        catch e
            @warn e "fail 3rd time, LAP_julia not loaded"
        end
    end
end

using ImageFiltering, Images, Interpolations, TestImages, ColorVectorSpace, FileIO

img = testimage("mandril_gray");

1
