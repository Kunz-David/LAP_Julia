
using ImageFiltering: centered, KernelFactors.gaussian, kernelfactors, imfilter
using LAP_julia
# using LAP_julia: Image

"""
    pad_images(image1::Image, image2::Image)

Adds zeros to the right and bottom of `image1` and `image2` to make them the same size.
"""
function pad_images(image1::Image, image2::Image)

    (a, b) = size(image1)
    (c, d) = size(image2)

    if (a < c)
        image1 = [image1; zeros(c - a, b)];
        a = c
    elseif (a > c)
        image2 = [image2; zeros(a - c, d)];
        c = a
    end

    if (b < d)
        image1 = [image1 zeros(a, d - b)];
    elseif (b > d)
        image2 = [image2 zeros(c, b - d)];
    end
    return image1, image2
end


"""
    rescale_intensities(image1::Image, image2::Image)

Rescale `image1` and `image2` intensities to span the whole `[0, 1]`.
"""
function rescale_intensities(image1::Image, image2::Image)

    max_intensity = maximum([image1[:]; image2[:]])
    min_intensity = minimum([image1[:]; image2[:]])

    image1 = (image1 .- min_intensity)./(max_intensity - min_intensity).*1;
    image2 = (image2 .- min_intensity)./(max_intensity - min_intensity).*1;

    return image1, image2
end

"""
    clean_using_gaussain(u::Matrix{<:Number}, window_half_size_one_dim::Integer)

Clean the Matrix `u` by smoothing using a square 2D Gaussian filter of size `2 * window_half_size_one_dim + 1` in each dimension.
"""
@inline function clean_using_gaussain(u::Matrix{<:Number}, window_half_size_one_dim::Integer)# where T<:Integer
    window_half_size = [window_half_size_one_dim, window_half_size_one_dim]
    return clean_using_gaussain(u, window_half_size)
end

"""
    clean_using_gaussain(u::Matrix{<:Number}, window_half_size)

Clean the Matrix `u` by smoothing using a 2D Gaussian filter of size `2 * window_half_size + 1`.
"""
function clean_using_gaussain(u::Matrix{<:Number}, window_half_size)

    σ_1 = 2 * window_half_size[1]
    σ_2 = 2 * window_half_size[2]

    cent_inds_1 = centered(-σ_1:σ_1)
    cent_inds_1 = centered(-σ_2:σ_2)

    gaus_1 = gaussian(σ_1, 2 * σ_1 + 1)
    gaus_2 = gaussian(σ_2, 2 * σ_2 + 1)

    kernf = kernelfactors((gaus_1, gaus_2))

    u_out = imfilter(u, kernf, "symmetric")
    return u_out
end

"""
    function mse(x, y)

Calculates the mean squared error between `x` and `y`.
"""
function mse(x, y)
    @assert (size(x) == size(y)) "sizes of $x and $y don't match"
    return mean((x .- y).^2)
end

"""
    function angle_rms(x, y)

Calculate the root mean square error in angle between `x` and `y`. Output in degrees.
"""
function angle_rms(x, y)
    @assert eltype(x) == eltype(y)
    @assert eltype(x) <: Complex

    return sqrt(mse(rad2deg.(angle.(x)), rad2deg.(angle.(y))))
end
"""
    function angle_mean(x, y)

Calculate the mean error in angle between `x` and `y`. Output in degrees.
"""
function angle_mean(x, y)
    return mean(abs.(rad2deg.(angle.(x)) - rad2deg.(angle.(y))))
end

"""
    function vec_len(x)

Calculate the lenght of vector `x`. `x` is a complex number.
"""
function vec_len(x)
    return sqrt(real(x)^2 + imag(x)^2)
end

"""
    function mean(x)

Calculate the mean of `x`.
"""
function mean(x)
    sum(x)/length(x)
end


# function maxim(a)
#     maximum(x->isnan(x) ? -Inf : x, a)
# end
#
# function minim(a)
#     minimum(x->isnan(x) ? +Inf : x, a)
# end
