module interpolation

export warp_img, interpolate_flow

using LAP_julia
using Interpolations, ColorVectorSpace, ImageFiltering, Images, ScatteredInterpolation

"""
    function warp_img(img, dx, dy; border_strat::Symbol=:replicate)

Warp the image `img` by `dx` in *x* direction and by `dy` in *y* direction.

`border_strat` indicates how to act when the resulting coordinate outside of bounds of `img`. Values: :replicate, :zeros.

See also: [`showflow`](@ref), [`imgshowflow`](@ref), [`Flow`](@ref)
"""
function warp_img(img, dx, dy; border_strat::Symbol=:replicate)
    itp = Interpolations.interpolate(img, BSpline(Linear()))
    if border_strat == :replicate
        etp = Interpolations.extrapolate(itp, Flat()) # added extrapolation
    elseif border_strat == :zeros
        etp = Interpolations.extrapolate(itp, 0)
    end
    inds = indices_spatial(img)
    rng = extrema.(inds)
    imw = similar(img, eltype(itp))
    for I in CartesianIndices(inds)
        y, x = Tuple(I) .+ (dy[I], dx[I])
        imw[I] = etp(y, x)
    end
    return imw
end

meshgrid(x, y) = [repeat(x, outer=length(y)) repeat(y, inner=length(x))]
meshgrid(x::Real, y::Real) = [repeat(1:x, outer=y) repeat(1:y, inner=x)]

"""
    interpolate_flow(flow::Flow, inds::Array{CartesianIndex, 1}, method::T=Multiquadratic(2)) where {T <: ScatteredInterpolation.RadialBasisFunction}

Interpolate a flow using the a Multiquadratic RBF with the parameter `ε` and using only the values at `inds` for the interpolation.

# Examples of `method` from ScatteredInterpolation:
Multiquadratic(ɛ = 1)
```math
ϕ(r) = \\sqrt{1 + (ɛr)^2}
```
Polyharmonic(k = 1)
```math
ϕ(r) = r^k, k = 1, 3, 5, ...
\\\\
ϕ(r) = r^k ln(r), k = 2, 4, 6, ...
```
(See ScatteredInterpolation for more.)

See also: [`showflow`](@ref), [`Flow`](@ref)
"""
function interpolate_flow(flow::Flow, inds::Array{CartesianIndex, 1}, method::T=Multiquadratic(2)) where {T <: ScatteredInterpolation.RadialBasisFunction}
    # get samples of flow at points
    samples = [flow[ind] for ind in inds]
    non_nan = filter(x -> !isnan(samples[x]), 1:length(samples))
    @assert length(non_nan) != 0

    # interpolate
    itp = ScatteredInterpolation.interpolate(method, LAP_julia.inds_to_points(inds[non_nan]), samples[non_nan]);
    gridPoints = meshgrid(size(flow)...)'
    interpolated = ScatteredInterpolation.evaluate(itp, gridPoints)

    gridded = reshape(interpolated, size(flow))

    return gridded
end

end #module
