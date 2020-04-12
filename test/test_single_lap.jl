
tile_size = 50
board_size = 4 # must be even
mini_board = [zeros(tile_size, tile_size) ones(tile_size, tile_size);
              ones(tile_size, tile_size) zeros(tile_size, tile_size)]

chessboard = repeat(mini_board, outer=(convert(Integer,(board_size/2)), convert(Integer,(board_size/2))))
img = chessboard

img[55, 40] = 0;
img[75, 40] = 0;
img[65, 30] = 0;
img[65, 20] = 0;

real_coef = 5;
imag_coef = -2;
flow = real_coef .* ones(size(img)) .+ imag_coef .* im .* ones(size(img))

img_showflow(img, flow)

imgw = LAP_julia.interpolation.imWarp_replicate(img, real(flow), imag(flow))
imgshow(imgw)

img_showflow(imgw, zeros(size(imgw)))
img_showflow(img, zeros(size(imgw)))
save("imgw.png", reverse(imgw, dims=1))
save("img.png", reverse(img, dims=1))


img_showflow(imgw, flow)

u_est, coeffs = LAP_julia.lap.single_lap(img, imgw, 3, 7, [15, 15])

showflow(u_est)
showflow(flow)
u_est[55,49]

## try to use opposite of u_est to get img from imgw


img_remade = LAP_julia.interpolation.imWarp(imgw, real(u_est), imag(u_est))

fig = PyPlot.figure(dpi = 300, figsize = (5, 5));
PyPlot.imshow(img_remade, cmap = :gray);
ax = gca();
ax.set_ylim(0, size(img_remade)[1]);
gcf()

save("img_remade.png", reverse(img_remade, dims=1))

fig = PyPlot.figure(dpi = 300, figsize = (5, 5));
PyPlot.imshow(img, cmap = :gray);
ax = gca();
ax.set_ylim(0, size(img)[1]);
gcf()

# difference:
img_dif = img - img_remade;
img_dif[40,101]

fig = PyPlot.figure(dpi = 300, figsize = (5, 5));
PyPlot.imshow(img_dif, cmap = :gray);
ax = gca();
ax.set_ylim(0, size(img_dif)[1]);
gcf()

## inpaint u_est

u_inp = u_est;
LAP_julia.inpaint.inpaint_nans!(u_inp)

showflow(u_inp)
showflow(flow)

mean(real(flow))
mean(imag(flow))

mean(real(u_inp))
mean(imag(u_inp))

mean(real(flow))
mean(imag(flow))

img_inp = LAP_julia.interpolation.imWarp(imgw, real(u_inp), imag(u_inp))

imgshow(img_inp)
imgshow(img)
imgshow(imgw)

save("img_inp.png", reverse(img_inp, dims=1))

# difference:
img_dif2 = img - img_inp;

fig = PyPlot.figure(dpi = 300, figsize = (5, 5));
PyPlot.imshow(img_dif2, cmap = :gray);
ax = gca();
ax.set_ylim(0, size(img_dif2)[1]);
gcf()

img_dif2[50, 75]
