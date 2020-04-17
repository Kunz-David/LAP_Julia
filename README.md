# LAP_julia

Local all-pass filtering registration method implemented in Julia from [this paper](http://www.ee.cuhk.edu.hk/~tblu/monsite/pdfs/gilliam1701.pdf).

## How to get use this module

Open up a Julia terminal and type:
```Julia
using Pkg; Pkg.add(PackageSpec(url="https://github.com/Kunz-David/LAP_Julia"))
```

----

## Algorithm functions
1) **`u_est`, `source_reg`, `figs` = `polyfilter_lap(target, source)`**
    - Implements the basic concept of `Algorithm 2` from the [paper](http://www.ee.cuhk.edu.hk/~tblu/monsite/pdfs/gilliam1701.pdf) without some features.
    - Finds `u_est` that transforms the image `source` image to the `target` image.
    - It uses `single_lap` iteratively; in each iteration using the transformation estimated by `single_lap` to warp the `source` image closer to the `target` image and then using this warped closer image as the source image in the next iteration, while using progressively smaller `filter_half_sizes` to estimate even small and faster varying displacements.
    - ### Input:
        - `target` is an image we want `source` to look like.
        - `source` is an warped image we want to transform into `target`.
    - ### Output:
        - `u_est` is the complex vector field that transforms `source` closer to `target`.
        - `source_reg` is the image `source` transformed by `u_est`.
        - `figs` is a 2D array of `PyPlot` Figures which shows the work of the algorithm at each iteration. For each iteration there are 3 Figures in this order: 1) current `u_est`, 2) newest addition to `u_est` `Δ_u`, 3) current `source_reg`.

2) **`u_est`, `all_coeffs` = `single_lap(image_1, image_2, filter_num::Int, filter_half_size::Int, window_size)`**
    - Implements the `Algorithm 1` from this [paper](http://www.ee.cuhk.edu.hk/~tblu/monsite/pdfs/gilliam1701.pdf).
    - Uses a gaussian filter basis to estimate a smoothly varying displacement field that transforms `image_2` to `image_1`.
    - ### Input:
        - `image_1` is a grayscale target image.
        - `image_2` is a grayscale source image.
        - `filter_num` number of basis filters used (so far only =3 implemented).
        - `filter_half_size` is the half size of the base of the gaussian filters used. (filter_size = 2*filter_half_size + 1) filter size: `filter_size` $\times$ `filter_size`.
        - `window_size` is the size of the local window (list of 2 ints) usually same as filter_size
    - ### Output:
        - `u_est` a complex vector field of size of `image_1` which is the displacement that transforms `image_2` closer to `image_1`.
        - `all_coeffs` are the coefficients of the basis filters for each pixel.

## Visualisation functions
1) **`imgshow(img)`**
    - Plots the image `img` putting the origin in the bottom left corner.
    - _Note:_ Use `img = reverse(img, dims=1)` on images from the `TestImages` module.
2) **`showflow(u_est, skip_count=nothing; fig=nothing, mag=1, legend=true)`**
    - Plots the complex flow `u_est`.
    - ### Input:
        - `u_est` is the complex flow.
        - _Optional:_ `skip_count` is the number of vectors to skip for the displaying of larger flows to be meaningful. Will be set to a magic constant depending on the size of `u_est`.
        - _Optional:_ `mag` is the magnification of the vectors.
3) **`img_showflow(imgw, flow; skip_count=nothing, mag=1)`**
    - A combination of `imgshow` and `showflow`.

## Other functions
1) **`inpaint_nans!(u)`**
    - Takes a complex matrix `u` containing _some_ NaNs and fills the NaNs using [this algorithm](https://www.researchgate.net/publication/220903053_Fast_Digital_Image_Inpainting)
2) **`imw` = `function imWarp(img, dx, dy)`**
    - Warps image `img` by `dx` in **x** direction and by `dy` in **y** direction and returns the warped image `imw`.
    - Sets pixels mapped out of the original image to 0.
3) **`imw` = `imWarp_replicate(img, dx, dy)`**
    - Warps image `img` by `dx` in **x** direction and by `dy` in **y** direction and returns the warped image `imw`.
    - Sets pixels mapped out of the original image to the closest pixel value still in `img`.
4) **`rand_flow` = `gen_rand_flow(flow_size, max_magnitude, tile_size=Nothing, filter_amp=Nothing)`**
    - Generates a smoothly varying flow `rand_flow` using uniform distribution.
    - ### Inputs:
        - `flow_size` is a tuple of dimensions of the flow.
        - `max_magnitude` upper bound to the amplitude of the displacement.
        - _Optional:_ `tile_size` size random element of start matrix (the larger the slower the flow). _Note_: If set larger than the `flow_size` it will generate a uniform pixel shift in a random direction.
        - _Optional:_ `filter_amp` is the size of the gaussian filter which is used to smooth the random start matrix
    - ### Output:
        - `rand_flow` is the generated random flow.


-------
## Basic usage example

1) **Make example images**

    chessboard
    ```Julia
    tile_size = 50
    board_size = 4 # must be even
    mini_board = [zeros(tile_size, tile_size) ones(tile_size, tile_size);
                  ones(tile_size, tile_size) zeros(tile_size, tile_size)]
    chessboard = repeat(mini_board, outer=(convert(Integer,(board_size/2)), convert(Integer,(board_size/2))))

    img = chessboard
    ```

    mandrill
    ```Julia
    img = testimage("mandril_gray")
    imgr = reverse(img, dims=1)

    img = imgr
    ```

2) **Generate a random flow**

    smoothly varying
    ```Julia
    flow = gen_rand_flow(size(img), 30);
    showflow(flow)
    ```

    uniform pixel shift
    ```Julia
    flow = gen_rand_flow(size(img), 30000);
    showflow(flow)
    ```

3) **Warp the input image**

    ```Julia
    imgw = LAP_julia.interpolation.imWarp_replicate(img, real(flow), imag(flow));
    imgshow(imgw)
    ```

4) **Test algorithm functions**

    single_lap
    ```Julia
    u_est, coeffs = single_lap(img, imgw, 3, 32, [65, 65]);
    LAP_julia.inpaint.inpaint_nans!(u_est);

    dewarped_img = LAP_julia.interpolation.imWarp_replicate(imgw, real(u_sin_est), imag(u_est))
    ```

    polyfilter_lap
    ```Julia
    u_est, source_reg, figs = LAP_julia.lap.polyfilter_lap(img, imgw);
    ```

5) **Show results**

    single_lap
    ```Julia
    imgshow(dewarped_img)
    showflow(u_est)
    ```

    polyfilter_lap
    ```Julia
    imgshow(source_reg)
    showflow(u_est)
    ```

    estimated flow comparison
    ```Julia
    showflow(flow) # original
    showflow(u_est .* (-1)) # estimated

    # difference
    showflow(u_est .* (-1) .- flow)
    ```
