module helpers

using ImageFiltering: centered, KernelFactors.gaussian, kernelfactors, imfilter

"""
    pad_images(image1, image2)

Adds zeros to the right and bottom of `image1` and `image2` to make them the same size.
"""
function pad_images(image1, image2)

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
    rescale_intensities(image1, image2)

Rescales `image1` and `image2` intensities to span the whole [0, 1].
"""
function rescale_intensities(image1, image2)

    max_intensity = maximum([image1[:]; image2[:]])
    min_intensity = minimum([image1[:]; image2[:]])

    image1 = (image1 .- min_intensity)./(max_intensity - min_intensity).*1;
    image2 = (image2 .- min_intensity)./(max_intensity - min_intensity).*1;

    return image1, image2
end

"""
    clean_using_gaussain(u, window_half_size)

function cleans the estimate of the displacement field `u` by smoothing using a Gaussian filter
"""
function clean_using_gaussain(u, window_half_size)

    # Sigma value for Gaussian filter:
    if length(window_half_size) == 1
        window_half_size = [window_half_size, window_half_size]
    end
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

end # module
