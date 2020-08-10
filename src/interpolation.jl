
using Interpolations, ColorVectorSpace, ImageFiltering, Images, ScatteredInterpolation


#TODO edit doc
"""
    warp_img(img::Image, dx, dy; border_strat::Symbol=:replicate)::Image
    warp_img(img::Image, dx, dy, filling_img::Image)::Image

Warp the image `img` by `dx` in *x* direction and by `dy` in *y* direction.
`border_strat` indicates how to act when the resulting coordinate outside of bounds of `img`.

# Possible variations:
- `border_strat == :replicate`: fill with the closest pixel value.
- `border_strat == :zeros`: fill with zeros.
- supply a `filling_img` as 3rd argument: fill with the provided image `filling_img`

See also: [`showflow`](@ref), [`imgshowflow`](@ref), [`Flow`](@ref)
"""
function warp_img(img, dx, dy; border_strat::Symbol=:replicate)
    if border_strat == :replicate
        return img_warp_replicate(img, dx, dy)
    elseif border_strat == :zeros
        return warp_img_zeros(img, dx, dy)
    end
end

function warp_img(img::Image, dx, dy, filling_img::Image)::Image
    return warp_img_with_filling(img, dx, dy, filling_img)
end


function warp_img_zeros(img::Image, dx, dy)::Image
    itp = Interpolations.interpolate(img, BSpline(Linear()))
    etp = Interpolations.extrapolate(itp, 0)
    inds = indices_spatial(img)
    rng = extrema.(inds)
    imgw = similar(img, eltype(itp))
    for I in CartesianIndices(inds)
        y, x = Tuple(I) .+ (dy[I], dx[I])
        imgw[I] = etp(y, x)
    end
    return imgw
end

function warp_img_with_filling(img::Image, dx, dy, filling_img::Image)::Image
    itp = Interpolations.interpolate(img, BSpline(Linear()))
    inds = indices_spatial(img)
    rng = extrema.(inds)
    imgw = similar(img, eltype(itp))
    for I in CartesianIndices(inds)
        imgw[I] = interplate_or_fill(I, dx[I], dy[I], rng, itp, filling_img)
    end
    return imgw
end

@inline function interplate_or_fill(ind, dxi, dyi, rng, itp, filling_img)
    y, x = (ind[1]+dyi, ind[2]+dxi)
    if is_in_bounds(y, rng[1]...) && is_in_bounds(x, rng[2]...)
        return itp(y, x)
    else
        return filling_img[ind]
    end
end

@inline function is_in_bounds(num, lo, hi)
    return lo <= num <= hi
end

function img_warp_replicate(img, dx, dy)
    itp = Interpolations.interpolate(img, BSpline(Linear()))
    inds = indices_spatial(img)
    rng = extrema.(inds)
    imw = similar(img, eltype(itp))
    for I in CartesianIndices(inds)
        dxi, dyi = dx[I], dy[I]
        y, x = clamp(I[1]+dyi, rng[1]...), clamp(I[2]+dxi, rng[2]...)
        imw[I] = itp(y, x)
    end
    return imw
end


meshgrid(x, y) = [repeat(x, outer=length(y)) repeat(y, inner=length(x))]
meshgrid(x::Real, y::Real) = [repeat(1:x, outer=y) repeat(1:y, inner=x)]

"""
interpolate_flow(flow_at_inds,
                 inds::Array{CartesianIndex{2}, 1},
                 flow_size;
                 method::Symbol=:quad,
                 kwargs=Dict())

Interpolate a complex dispalcment field of size `flow_size` from an array of displacement vectors `flow_at_inds`, that are at the locations `inds`.
You can choose between two methods `method=:quad` which fits the displacement vectors into a global quadratic displacement field
and `method=:rbf` which uses an rbf model. Additional keyword arguments of the chosen method can be passed as `kwargs`.

See also: [`Flow`](@ref), [`interpolate_flow_quad`](@ref), [`interpolate_flow_rbf`](@ref)
"""
function interpolate_flow(flow_at_inds,
                          inds::Array{CartesianIndex{2}, 1},
                          flow_size;
                          method::Symbol=:quad,
                          kwargs=Dict())

if method == :quad
    generated_flow = interpolate_flow_quad(flow_at_inds, inds, flow_size)
elseif method == :rbf
    generated_flow = interpolate_flow_rbf(flow_at_inds, inds, flow_size; kwargs...)
end
    return generated_flow
end


#TODO edit docs
"""
    interpolate_flow_rbf(flow_at_inds,
                         inds::Array{CartesianIndex{2}, 1},
                         flow_size,
                         rbf_method::T=Multiquadratic(2)) where {T <: ScatteredInterpolation.RadialBasisFunction}

Interpolate a complex displacement field of size `flow_size` from displacement vectors `flow_at_inds` at `inds`
using the a Multiquadratic RBF model with the parameter `ε` and using only the values at `inds` for the interpolation.

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

See also: [`showflow`](@ref), [`Flow`](@ref), [`interpolate_flow`](@ref), [`interpolate_flow_quad`](@ref)
"""
function interpolate_flow_rbf(flow_at_inds,
                              inds::Array{CartesianIndex{2}, 1},
                              flow_size,
                              rbf_method::T=Multiquadratic(2)) where {T <: ScatteredInterpolation.RadialBasisFunction}

    non_nan_locs = filter(x -> !isnan(flow_at_inds[x]), 1:length(flow_at_inds))
    @assert length(non_nan_locs) != 0

    # interpolate
    itp = ScatteredInterpolation.interpolate(rbf_method, inds_to_points(inds[non_nan_locs]), flow_at_inds[non_nan_locs]);
    gridPoints = meshgrid(flow_size...)'
    interpolated = ScatteredInterpolation.evaluate(itp, gridPoints)

    gridded = reshape(interpolated, flow_size)

    return gridded
end

"""
    interpolate_flow_quad(flow_at_inds, inds, flow_size)

Interpolate a complex displacement field of size `flow_size` from displacement vectors `flow_at_inds` at `inds`
using a global quadratic model.

The fit to this model is made by minimizing the following:

```math
\\min_{a_{k}} \\sum_{x, y ∈ \\mathrm{inds}}\\left|u_{\\mathrm{quad}}(x, y)-u_{\\mathrm{flow\\_at\\_inds}}(x, y)\\right|^{2}
```,

where:

```math
u_{\\mathrm{quad}}(x, y) = a_1 + a_2x + a_3y + a_4x^2 + a_5y^2 + a_6xy
```


See also: [`showflow`](@ref), [`Flow`](@ref), [`interpolate_flow`](@ref), [`interpolate_flow_rbf`](@ref)
"""
function interpolate_flow_quad(flow_at_inds, inds, flow_size)
    A = make_A(inds)
    coeffs = qr(A, Val(true))\flow_at_inds
    estim_flow = make_flow2(coeffs, flow_size)
end

function make_flow2(coeffs, flow_size)
    # flow = similar(coeffs, ComplexF64, flow_size)
    flow = Array{ComplexF64}(undef, flow_size)
    xs = collect(1:flow_size[1])
    if flow_size[1] == flow_size[2]
        ys = xs
    else
        ys = collect(1:flow_size[2])
    end
    flow .= coeffs[1]
    flow .+= (xs .* coeffs[2])
    transpose(flow) .+= (ys .* coeffs[3])
    flow .+= ((xs.^2) .* coeffs[4])
    transpose(flow) .+= ((ys.^2) .* coeffs[5])
    flow .+= ((repeat(xs, 1, flow_size[2]) .* ys') .* coeffs[6])
    return flow
end

function make_A(inds)
    points = transpose(inds_to_points(inds))
    A = hcat(ones(length(inds)), points[:, 1], points[:,2], points[:,1].^2, points[:,2].^2, points[:,1] .* points[:,2])
end
