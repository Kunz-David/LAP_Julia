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
