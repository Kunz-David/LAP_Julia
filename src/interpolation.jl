module interpolation

using Interpolations, ColorVectorSpace, ImageFiltering, Images

"""
    function imWarp(img, dx, dy)

Warps image `img` by `dx` in *x* direction and by `dy` in *y* direction

- sets pixels mapped out of the original image to 0
- returns warped image
"""
function imWarp(img, dx, dy)
    itp = interpolate(img, BSpline(Linear()))
    inds = indices_spatial(img)
    rng = extrema.(inds)
    imw = similar(img, eltype(itp))
    for I in CartesianIndices(inds)
        dxi, dyi = dx[I], dy[I]
        y, x = I[1]+dyi, I[2]+dxi
        if is_in(y, rng[1]) && is_in(x,rng[2])
            imw[I] = itp(y, x)
        else
            imw[I] = 0
        end
    end
    return imw
end

function is_in(x, (low, high))
    return x >= low && x <= high
end

function imWarp_replicate(img, dx, dy)
    itp = interpolate(img, BSpline(Linear()))
    inds = indices_spatial(img)
    rng = extrema.(inds)
    imw = similar(img, eltype(itp))
    for I in CartesianIndices(inds)
        dxi, dyi = dx[I], dy[I]
        y, x = clamp(I[1]+dyi, rng[1]...), clamp(I[2]+dxi, rng[2]...)
        @assert (!isnan(y) && !isnan(x)) "are nans"
        @assert is_in(y, (rng[1])) "y is %i" y
        @assert is_in(x, (rng[2])) "x is %i" x
        imw[I] = itp(y, x)
    end
    return imw
end

end #module
