
using Images, ImageFiltering, LinearAlgebra


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

"""
    find_edge_points(img::Image; spacing::Int=10, number::Int=100, sigma=1, mask::BitArray{2}=[])

Greadily locate `number` of points in `img` with the highest gradient magnitude. These points have will be atleast `spacing`
pixels appart, they are returned as an array of `CartesianIndex`. A `mask` can be used to exclude some parts of `img` from the search.

Note: `sigma` is the gaussian filter used to smooth `img` before looking for gradients.
"""
function find_edge_points(img::Image;
                          spacing::Int=40,
                          number::Int=35,
                          sigma=1,
                          mask::BitArray{2}=BitArray{2}(undef, 0,0))

    # calculate gradient mag
    grad, mag = gradient_magnitude(ImageFiltering.imfilter(img, ImageFiltering.Kernel.gaussian(sigma)))

    @assert typeof(mag) <: Array

    #const
    marker = zero(eltype(mag)) # marker for masked gradients we don't want keypoints at the points of 0 gradient
    # marker = 0.25
    if !isempty(mask)
        @assert( size(mag) == size(mask))
        for i = 1:length(mag)
            if mask[i] == 0
                mag[i] = marker
            end
        end
    end

    N = 2
    ind = sortperm(reshape(mag,length(mag)), alg=PartialQuickSort(2*number*spacing^(N-1)),rev=true) # check here
    numkpts = 0 # number of keypoints found so far

    kpts = Array{CartesianIndex}(undef, number)

    # go through indices from the largest gradient
    for i in ind

        if numkpts >= number break; end
        if mag[i] == marker continue; end

        # @assert(mag[i] > 10 * eps(Float64))

        # good keypoint
        numkpts += 1

        pos = CartesianIndices(size(mag))[i]
        ppos = Tuple(pos) # position in pixel units

        kpts[numkpts] = pos

        # prevent the neighbors from being added
        fill_circle!(mag, ppos, spacing, marker)
    end # for

    if numkpts < number # shrink array if less keypoints were found
        resize!(kpts, numkpts)
    end
    return kpts
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
#         fill_circle!(mag, ppos, spacing, marker)
#     end # for i
#
#     return kpts, mag_save, mag
# end # function find_keypoints_from_gradients

"""
    fill_circle!(mag, pos, spacing, marker)

Fill a circle with a radius of `spacing` around the point `pos` into the matrix `mag` with `marker`.
"""
function fill_circle!(mag, pos, spacing, marker)
    mag_size = size(mag)
    for k = (pos[1]-spacing <= 1 ? 1 : pos[1]-spacing) : (pos[1]-spacing >= mag_size[1] ? pos[1]-spacing : mag_size[1]),
        l = (pos[2]-spacing <= 1 ? 1 : pos[2]-spacing) : (pos[2]-spacing >= mag_size[2] ? pos[2]-spacing : mag_size[2])
        # @assert (k >= 1) k
        if (sum(((k, l) .- pos).^2) <= spacing^2)
            mag[k, l] = marker
        end
    end
end
