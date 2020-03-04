
for k = 1:3, l = 1:3, m = 1:3, n = 1:3

        im1 = ones(k, l)
        im2 = ones(m, n)

        out = im1, im2
        size((LAP_julia.pad_images(im1, im2))[1])
end
