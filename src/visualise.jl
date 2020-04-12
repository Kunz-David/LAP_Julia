module visualise

using PyPlot: quiver, gcf, figure
using Printf

function imgshow(img)
    fig = PyPlot.figure(dpi = 300, figsize = (5, 5));
    PyPlot.imshow(img, cmap = :gray);
    ax = gca();
    ax.set_ylim(0, size(img)[1]);
    gcf()
end

"""
    showflow(u_est, skip_count=0)

shows displacement field `u_est` by default not skipping any vector.
The second optional parameter: `skip_count` determines how many vectors to skip in the visualisation.
"""
function showflow(u_est, skip_count=nothing; fig=nothing, mag=1, legend=true)

    # set defaults
    if skip_count == nothing
        magic_const = 20
        skip_count = ceil(Int64, maximum(size(u_est))/magic_const)-1
    end

    # prepare data
    uv_flow = zeros(size(u_est)..., 2)
    uv_flow[:, :, 1] = real(u_est)
    uv_flow[:, :, 2] = imag(u_est)

    n = skip_count+1

    siz = size(u_est)

    # meshgrid(x, y) = (repeat(x, outer=length(y)), repeat(y, inner=length(x)))

    # The x coordinates of the arrow locations
    X = [1:siz[2];]

    # The y coordinates of the arrow locations
    Y = [1:siz[1];]

    trimmed_real = fill(NaN, size(uv_flow[:,:,1]))
    trimmed_real[1:n:end, 1:n:end] .= uv_flow[1:n:end, 1:n:end, 1]
    trimmed_imag = fill(NaN, size(uv_flow[:,:,2]))
    trimmed_imag[1:n:end, 1:n:end] .= uv_flow[1:n:end, 1:n:end, 2]

    max = maximum(filter(!isnan, abs.([trimmed_imag trimmed_real])))

    scaling = max/(n * mag)
    color = atan.(trimmed_real, trimmed_imag)
    color = (color .- minimum(color)) ./ maximum(color)

    if fig == nothing
        fig, ax = subplots(dpi = 300)
    else
        fig = fig
        ax = gca()
    end

    q = ax.quiver(X, Y, trimmed_real, trimmed_imag, atan.(trimmed_real, trimmed_imag),
            scale_units = "xy", scale = scaling)

    # angles='uv' sets the angle of vector by atan2(u,v), angles='xy' draws the vector from (x,y) to (x+u, y+v)
    ax.set_aspect(1.)
    xlabel("x - real")
    ylabel("y - imag")
    title("Flow")

    if legend == true
        # subplots_adjust(bottom=0.1)
        label_text = "length = " * @sprintf("%0.2f", max)
        ax.quiverkey(q, X=0.19, Y = 0.04, U = max, coordinates="figure",
        label=label_text,labelpos = "E")
    end

    return gcf()
end

using PyPlot

function img_showflow(imgw, flow; skip_count=nothing, mag=1)

    fig = PyPlot.figure(dpi = 300, figsize = (5, 5));
    PyPlot.imshow(imgw, cmap = :gray);
    ax = gca();
    ax.set_ylim(0, size(imgw)[1]);

    showflow(flow, skip_count, fig=fig, mag=mag);
    gcf()

end

function warp_img_showflow(img, flow; skip_count=nothing, mag=1)
    imgw = LAP_julia.interpolation.imWarp(img, real(flow), imag(flow))

    img_showflow(imgw, flow, skip_count=skip_count, mag=mag)
end

end #module
