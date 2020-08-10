
# function better_overlay(img1, img2)
using PyPlot


a, b, c = rand(5, 5), rand(5, 5), rand(5, 5)


function img_overlay_v2(orange_img, blue_img)

    normed_orange_img = orange_img.*(1/extrema(orange_img)[2])
    normed_blue_img = blue_img.*(1/extrema(blue_img)[2])

    color1 = [1 0.7 0.5]
    color2 = [0 0.3 0.5]

    a = normed_orange_img
    b = color1[2].*normed_orange_img .+ color2[2].*normed_blue_img
    c = color1[3].*normed_orange_img .+ color2[3].*normed_blue_img

    imm = collect(colorview(RGB, a, b, c))
    return imm
end

normedimg = img.*(1/extrema(img)[2])
normedimgw = imgw.*(1/extrema(imgw)[2])

img
imgw

color1 = [0 0.7 0.5]
color2 = [0 0.3 0.5]

a = normedimg
b = color1[2].*normedimg .+ color2[2].*normedimgw
c = color1[3].*normedimg .+ color2[3].*normedimgw

imm = collect(colorview(RGB, a, b, c))
