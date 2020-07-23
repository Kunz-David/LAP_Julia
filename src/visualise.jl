
using PyPlot #: quiver, gcf, figure, imshow, subplots, title, xlabel, ylabel, gca
using Printf

# set plot size defaults
SMALL = 8
MEDIUM = 10
LARGE = 12

rc("font", size=SMALL)          # controls default text sizes
rc("axes", titlesize=LARGE)     # fontsize of the axes title
rc("axes", labelsize=SMALL)     # fontsize of the x and y labels
rc("xtick", labelsize=SMALL)    # fontsize of the tick labels
rc("ytick", labelsize=SMALL)    # fontsize of the tick labels
rc("legend", fontsize=SMALL)    # legend fontsize
rc("figure", titlesize=LARGE)   # fontsize of the figure title
rc("figure", figsize=(6,6))
rc("figure", dpi=400)
rc("savefig", dpi=400)
rc("image", origin="lower")
# rcdefaults()

# TODO: add docs
function addpoints(inds; ret::Symbol=:figure, labels=[])
    pos_x = [ind[1] for ind in inds]
    pos_y = [ind[2] for ind in inds]

    ax = PyPlot.scatter(pos_y, pos_x, marker = :x)

    for (k, label) in enumerate(labels)
        annotate(label, (pos_y[k], pos_x[k]))
    end

    if ret == :figure
        return gcf()
    elseif ret == :pyobject
        return ax
    end
end

# TODO: add docs
function imgoverlay(img1, img2; fig=nothing, figtitle::String="Image overlay", ret::Symbol=:figure)

    if fig == nothing
        fig, ax = subplots(dpi = 400)
    else
        fig = fig
        ax = gca()
    end

    img1 = Float64.(img1)
    img2 = Float64.(img2)
    imshow(img1, cmap = :Blues_r, alpha = 0.5)
    imshow(img2, cmap = :Oranges_r, alpha = 0.5)


    title(figtitle)

    if ret == :figure
        return gcf()
    elseif ret == :pyobject
        return ax
    end
end


"""
    imgshow(img; <keyword arguments>)

Return a figure with image `img`. Plot origin is in the botom left.

# Keyword Arguments
- `fig=nothing`: add a figure to plot in. By defaults creates a blank new figure.
- `figtitle::String="Image"`: add title to the figure.
- `ret::Symbol=:figure`: set return object, by default returns Figure, other options: :pyobject returns a PyObject. (Using figure makes Juno directly plot.)
- `origin_left_bot::Bool=false`: set the origin to the lower left corner, which is not typical for images (they will appear flipped).

See also: [`showflow`](@ref), [`imgshowflow`](@ref), [`warp_imgshowflow`](@ref)
"""
function imgshow(img;
                 fig=nothing,
                 figtitle::String="Image",
                 ret::Symbol=:figure,
                 origin_left_bot::Bool=false)

    if fig == nothing
        fig, ax = subplots(dpi = 400)
    else
        fig = fig
        ax = gca()
    end

    img = Float64.(img)

    imshow(img, cmap = :gray);
    # ax.set_ylim(0, size(img)[1]-1);
    # ax.set_xlim(0, size(img)[2]-1);
    if origin_left_bot == true
        ax.invert_yaxis()
    end

    title(figtitle)

    if ret == :figure
        return gcf()
    elseif ret == :pyobject
        return ax
    end
end

"""
    showflow(flow::Flow; <keyword arguments>)

Return a figure with the displacement field `flow` by default skipping some vectors to make it easy to read.

# Keyword Arguments
- `flow::Flow`: the vector flow to be plotted.
- `disp_type::Symbol=:full` : display mode, either `:full` -> display all the vectors, or `:sparse` -> displays only the vectors above a threshold, or `:auto` which decide for you based on the data.
- `skip_count=nothing`: the number of vectors to skip between each displayed vector. By default set so that the output is ``20 × 20`` vectors.
- `fig=nothing`: add a figure to plot in. By defaults creates a blank new figure.
- `mag::Real=1`: magnify the plotted vectors.
- `key::Bool=true`: add key with maximum vector length.
- `figtitle::String="Flow"`: add title to the figure.
- `ret::Symbol=:figure`: set return object, by default returns Figure, other options: :pyobject returns a PyObject. (Using figure makes Juno directly plot.)

See also: [`imgshow`](@ref), [`imgshowflow`](@ref), [`warp_imgshowflow`](@ref)
"""
function showflow(flow::Flow; disp_type::Symbol=:auto, skip_count=nothing, fig=nothing, mag::Real=1, key::Bool=true, figtitle::String="Flow", ret::Symbol=:figure)

    # set defaults
    vecs_in_one_dim = 20
    if skip_count == nothing
        skip_count = ceil(Int64, maximum(size(flow))/vecs_in_one_dim)-1
    end

    # prepare data
    uv_flow = zeros(size(flow)..., 2)
    uv_flow[:, :, 1] = real(flow)
    uv_flow[:, :, 2] = imag(flow)

    n = skip_count+1
    # if n == 1
    #     n = 2
    # end


    # meshgrid(x, y) = (repeat(x, outer=length(y)), repeat(y, inner=length(x)))
    # The x coordinates of the arrow locations
    X = [1:size(flow, 2);]

    # The y coordinates of the arrow locations
    Y = [1:size(flow, 1);]

    min_threshold = 1e-10

    if disp_type == :auto
        if count(x -> vec_len(x) > min_threshold, flow)/length(flow) <= 0.1
            disp_type = :sparse
        else
            disp_type = :full
        end
    end

    if disp_type == :full
        trimmed_real = fill(NaN, size(uv_flow[:,:,1]))
        trimmed_real[1:n:end, 1:n:end] .= uv_flow[1:n:end, 1:n:end, 1]
        trimmed_imag = fill(NaN, size(uv_flow[:,:,2]))
        trimmed_imag[1:n:end, 1:n:end] .= uv_flow[1:n:end, 1:n:end, 2]
    elseif disp_type == :sparse
        cond = ((abs.(real(flow)) .> min_threshold) .| (abs.(imag(flow)) .> min_threshold))
        trimmed_real = fill(NaN, size(uv_flow[:,:,1]))
        trimmed_real[cond] .= uv_flow[:, :, 1][cond]
        trimmed_imag = fill(NaN, size(uv_flow[:,:,2]))
        trimmed_imag[cond] .= uv_flow[:, :, 2][cond]
        if all(isnan, trimmed_real) || all(isnan, trimmed_imag)
            trimmed_imag = zeros(size(trimmed_imag))
            trimmed_real = zeros(size(trimmed_real))
        end
    end

    maxi = maximum(filter(!isnan, abs.([trimmed_imag trimmed_real])))

    scaling = maxi/(n * mag)
    # color = atan.(trimmed_real, trimmed_imag)
    # color = (color .- minimum(color)) ./ maximum(color)

    if fig == nothing
        fig, ax = subplots(dpi = 400)
    else
        fig = fig
        ax = gca()
    end

    q = ax.quiver(X, Y, trimmed_real, trimmed_imag, atan.(trimmed_real, trimmed_imag),
            scale_units = "xy", scale = scaling)

    # angles='uv' sets the angle of vector by atan2(u,v), angles='xy' draws the vector from (x,y) to (x+u, y+v)
    ax.set_aspect(1.)
    xlabel("x - real\n\n")
    ylabel("y - imag")
    title(figtitle)

    # ax.invert_yaxis()

    if key == true
        # subplots_adjust(bottom=0.1)
        max_mag = max_displacement(flow)
        label_text = "length = " * @sprintf("%0.2f", max_mag)
        ax.quiverkey(q, X=0.75, Y = 0.2, U = max_mag, coordinates="inches",
        label=label_text, labelpos = "E")
    end

    if ret == :figure
        return gcf()
    elseif ret == :pyobject
        return q
    end
end


"""
    imgshowflow(img, flow; <keyword arguments>)

Return a figure with an image `img` and displacement field `flow`.

# Keyword Arguments
- `img`: the image to be plotted.
- `flow::Flow`: the vector flow to be plotted.
- `disp_type::Symbol=:full` : display mode, either `:full` -> display all the vectors, or `:sparse` -> displays only the vectors above a threshold.
- `skip_count=nothing`: the number of vectors to skip between each displayed vector. By default set so that the output is ``20 × 20`` vectors.
- `fig=nothing`: add a figure to plot in. By defaults creates a blank new figure.
- `mag::Real=1`: magnify the plotted vectors.
- `key::Bool=true`: add key with maximum vector length.
- `figtitle::String="Iamge with Flow"`: add title to the figure.
- `ret::Symbol=:figure`: set return object, by default returns Figure, other options: :pyobject returns a PyObject. (Using figure makes Juno directly plot.)
- `origin_left_bot::Bool=true`: set the origin to the lower left corner, which is not typical for images (they will appear flipped).

See also: [`imgshow`](@ref), [`showflow`](@ref), [`warp_imgshowflow`](@ref)

"""
function imgshowflow(img,
                     flow::Flow;
                     disp_type::Symbol=:full,
                     skip_count=nothing, fig=nothing,
                     mag::Real=1,
                     key::Bool=true,
                     figtitle::String="Image with Flow",
                     ret::Symbol=:figure,
                     origin_left_bot::Bool=true)

    fig = imgshow(img, fig=fig, figtitle="", ret=:figure, origin_left_bot=origin_left_bot)
    showflow(flow, disp_type=disp_type, skip_count=skip_count, fig=fig, mag=mag, key=key, figtitle=figtitle, ret=ret);
end

"""
    warp_imgshowflow(img, flow; <keyword arguments>)

Return a figure with an image and a displacement field, where the image is warped by the displacement field.

# Keyword Arguments
- `img`: the image to be plotted.
- `flow::Flow`: the vector flow to be plotted.
- `disp_type::Symbol=:full` : display mode, either `:full` -> display all the vectors, or `:sparse` -> displays only the vectors above a threshold.
- `skip_count=nothing`: the number of vectors to skip between each displayed vector. By default set so that the output is ``20 × 20`` vectors.
- `fig=nothing`: add a figure to plot in. By defaults creates a blank new figure.
- `mag::Real=1`: magnify the plotted vectors.
- `key::Bool=true`: add key with maximum vector length.
- `figtitle::String="Iamge with Flow"`: add title to the figure.
- `ret::Symbol=:figure`: set return object, by default returns Figure, other options: :pyobject returns a PyObject. (Using figure makes Juno directly plot.)
- `origin_left_bot::Bool=true`: set the origin to the lower left corner, which is not typical for images (they will appear flipped).

See also: [`imgshow`](@ref), [`showflow`](@ref), [`imgshowflow`](@ref)
"""
function warp_imgshowflow(img,
                          flow::Flow;
                          disp_type::Symbol=:full,
                          skip_count=nothing,
                          fig=nothing, mag::Real=1,
                          key::Bool=true,
                          figtitle::String="Warped Image with Flow",
                          ret::Symbol=:figure,
                          origin_left_bot::Bool=true)

    imgw = warp_img(img, -real(flow), -imag(flow))
    imgshowflow(imgw, flow, disp_type=disp_type, skip_count=skip_count, fig=fig, mag=mag, key=key, figtitle=figtitle, ret=ret, origin_left_bot=origin_left_bot);
end
