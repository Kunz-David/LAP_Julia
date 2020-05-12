module interpolation

export warp_img
using Interpolations, ColorVectorSpace, ImageFiltering, Images, LAP_julia

"""
    function warp_img(img, dx, dy; border_strat::Symbol=:replicate)

Warp the image `img` by `dx` in *x* direction and by `dy` in *y* direction.

`border_strat` indicates how to act when the resulting coordinate outside of bounds of `img`. Values: :replicate, :zeros.

See also: [`showflow`](@ref), [`imgshowflow`](@ref), [`Flow`](@ref)
"""
function warp_img(img, dx, dy; border_strat::Symbol=:replicate)
    itp = interpolate(img, BSpline(Linear()))
    if border_strat == :replicate
        etp = extrapolate(itp, Flat()) # added extrapolation
    elseif border_strat == :zeros
        etp = extrapolate(itp, 0)
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

end #module
