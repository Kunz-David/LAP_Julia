
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

# TODO add docs
function interpolate_flow(flow_at_inds,
                          inds::Array{CartesianIndex, 1},
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
    interpolate_flow_rbf(flow::Flow, inds::Array{CartesianIndex, 1}, method::T=Multiquadratic(2)) where {T <: ScatteredInterpolation.RadialBasisFunction}

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
function interpolate_flow_rbf(flow_at_inds,
                              inds::Array{CartesianIndex, 1},
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


function interpolate_flow_quad(flow_at_inds, inds, flow_size)
    A = make_A(inds)
    coeffs = A\flow_at_inds
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
    flow .+= ((xs' .* ys) .* coeffs[6])
    return flow
end

function make_A(inds)
    points = transpose(inds_to_points(inds))
    A = hcat(ones(length(inds)), points[:, 1], points[:,2], points[:,1].^2, points[:,2].^2, points[:,1] .* points[:,2])
end
