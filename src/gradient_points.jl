module gradient_points

export find_keypoints_from_gradients

using StatsBase
using Images, ImageFiltering, LinearAlgebra
using LAP_julia

"""
    gradient_magnitude(f)
Return the directional derivatives (components of the gradient)
and the gradient mag for both color and grayscale images. For
`spatialorder(f)=="xy"`, `grad[1]` is a derivative wrt x, the first
coordinate and `grad[2]` wrt to the second...

"""
# function gradient_magnitude(f::Image{T,2}) where {T<:ColorTypes.Gray}
#     @assert(length(size(f))==2) # only 2D supported for the moment
#     s=ImageFiltering.Kernel.sobel()     # vytvorime Sobelovsky filtr ve smeru x a y
#     grad=[Images.imfilter(f,s[i]) for i in 1:2]   # filtrujeme ve smeru x a y
#     return grad,Images.mag(grad[1],grad[2]) # vracime gradient a jeho normu
# end
"""
    gradient_magnitude(f)
Return the directional derivatives (components of the gradient)
and the gradient mag for both color and grayscale images. For
`spatialorder(f)=="xy"`, `grad[1]` is a derivative wrt x, the first
coordinate and `grad[2]` wrt to the second...

"""
function gradient_magnitude(f)
    @assert(length(size(f))==2) # only 2D supported for the moment
    s=ImageFiltering.Kernel.sobel()     # vytvorime Sobelovsky filtr ve smeru x a y
    grad=[Images.imfilter(f,s[i]) for i in 1:2]   # filtrujeme ve smeru x a y
    return grad, Images.magnitude(grad[1],grad[2]) # vracime gradient a jeho normu
end


# V dalši funkci hledáme tzv. klíčové body, tedy body ve kterých budeme později hledat neznámou deformaci (posunutí) ve vhodném okénku (ROI)


"""
    find_keypoints_from_gradients{T,N}(f::Image{T,N}; spacing = 10, numbe = 100, sigma = ones(Float64, N), nlength = 1., mask = [])
Take a color or grayscale image and return `SimpleKeypoint`s located at high-gradient positions
of the image `f`. Returns `number` keypoints at least `spacing` pixels apart. If this is not
possible, returns less keypoints. `sigma` is the Gaussian filter size. You can suppress keypoint
creation at some locations by setting particular voxels in `mask` to 0. Sizes of `mask` and `f`
should be equal.

The coordinates are in physical units.

"""
function find_keypoints_from_gradients(f::Image; spacing=10, number::Int=100, sigma=ones(Float64,2), nlength=1., mask=[])

    # calculate gradient mag
    grad, mag = gradient_magnitude(ImageFiltering.imfilter(f, ImageFiltering.Kernel.gaussian(sigma)))

    mag_save = copy(mag)

    @assert typeof(mag) <: Array

    #const
    marker = zero(eltype(mag)) # marker for masked gradients we don't want keypoints at the points of 0 gradient
    # marker = 0.25
    if mask != []
        @assert( size(mag) == size(mask))
        for i = 1:length(mag)
            if mask[i] == 0
                mag[i] = marker
            end
        end
    end
    # sort indices by decreasing mag
    #ind=@debugtime( collect(eachindex(mag)), "Indices")
    #@debugtime( sort!(ind,by=i->mag[i],rev=true), "Sort")
    N = 2
    ind = sortperm(reshape(mag,length(mag)),alg=PartialQuickSort(2*number*spacing^(N-1)),rev=true)
    numkpts = 0 # number of keypoints found so far

    kpts = Array{SimpleKeypoint}(undef, number)
    pspacing = map(float, Images.pixelspacing(f))
    rspacing = spacing ./ pspacing
    # go through indices from the largest gradient and put each
    for i in ind
        # if we already have enough keypoints
        if numkpts >= number break;    end
        # if this point is too close to another
        if mag[i] == marker continue
        end
        @assert(mag[i] > 10 * eps(Float64))
        # good keypoint, let"s add it
        numkpts += 1

        ppos = Tuple(CartesianIndices(size(mag))[i]) # position in pixel units
        pos = ([ppos...] .- 0.5) .* pspacing    # position in physical units

        norm = LinearAlgebra.norm(Float64[ grad[k][i] for k=1:N], 2) * nlength
        @assert (mag[i] ≈ norm) (mag[i], norm)
        kpts[numkpts] = SimpleKeypoint(ppos, norm)
        # prevent the neighbors from being added
        fill_ellipse!(mag, ppos, spacing, marker)
    end # for i

    if numkpts<number # shrink array if less keypoints were found
        resize!(kpts, numkpts)
    end
    return kpts, mag_save, mag
end # function find_keypoints_from_gradients





# function find_keypoints_from_gradients_p_field(f::Array{T, 2}; spacing=25, number::Int=100, sigma=ones(Float64,2), nlength=1., mask=[]) where {T}
#
#     # calculate gradient mag
#     grad, mag = gradient_magnitude(ImageFiltering.imfilter(f, ImageFiltering.Kernel.gaussian(sigma)))
#
#     mag_save = copy(mag)
#     mag = mag.^4
#
#     @assert typeof(mag) <: Array
#
#     #const
#     marker = zero(eltype(mag)) # marker for masked gradients we don't want keypoints at the points of 0 gradient
#     if mask != []
#         @assert( size(mag) == size(mask))
#         for i = 1:length(mag)
#             if mask[i] == 0
#                 mag[i] = marker
#             end
#         end
#     end
#     # sort indices by decreasing mag
#
#     kpts = Array{SimpleKeypoint}(undef, number)
#     # go through indices from the largest gradient and put each
#     for k in 1:number
#
#         ind = wsample(Array(1:length(mag)), mag[:])
#
#         # @assert(mag[ind] > 10 * eps(Float64))
#         # good keypoint, let"s add it
#
#         ppos = Tuple(CartesianIndices(size(mag))[ind]) # position in pixel units
#
#         N = 2
#         norm = LinearAlgebra.norm(Float64[ grad[l][ind] for l=1:N], 2) * nlength
#         @assert mag_save[ind] ≈ norm
#         kpts[k] = SimpleKeypoint(ppos, norm)
#         # prevent the neighbors from being added
#         fill_ellipse!(mag, ppos, spacing, marker)
#     end # for i
#
#     return kpts, mag_save, mag
# end # function find_keypoints_from_gradients
#

function fill_ellipse!(mag, pos, spacing, marker)
    mag_size = size(mag)
    for k = (pos[1]-spacing <= 1 ? 1 : pos[1]-spacing) : (pos[1]-spacing >= mag_size[1] ? pos[1]-spacing : mag_size[1]),
        l = (pos[2]-spacing <= 1 ? 1 : pos[2]-spacing) : (pos[2]-spacing >= mag_size[2] ? pos[2]-spacing : mag_size[2])
        @assert (k >= 1) k
        if (sum(((k, l) .- pos).^2) <= spacing^2)
            mag[k, l] = marker
        end
    end
end

end # module
