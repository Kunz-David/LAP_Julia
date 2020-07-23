
using FileIO, ImageIO, PyPlot

## lena reg intro
img, imgw, flow = gen_init()

imgshow(img, figtitle="Target")
savefig("../plots/intro_lena_orig.png")

imgshow(imgw, figtitle="Source")
savefig("../plots/intro_lena_warped.png")

imgoverlay(img, imgw, figtitle="Blending of Target and Source")
savefig("../plots/intro_lena_bl_or_movement.png")

showflow(flow, figtitle="Displacement field")
savefig("../plots/intro_lena_flow.png")


## lena reg intro with red rectangle
using PyCall
@pyimport matplotlib.patches as patches

function add_rect(ax)
    rect = patches.Rectangle((100,110),30,30,linewidth=1,edgecolor="r",facecolor="none")
    ax.add_patch(rect)
    return gcf()
end

img, imgw, flow = gen_init()

imgshow(img, figtitle="Target")
add_rect(gca())
savefig("../plots/rect_intro_lena_orig.png")

imgshow(imgw, figtitle="Source")
add_rect(gca())
savefig("../plots/rect_intro_lena_warped.png")

imgoverlay(img, imgw, figtitle="Blending of Target and Source")
add_rect(gca())
savefig("../plots/rect_intro_lena_bl_or_movement.png")

showflow(flow, figtitle="Displacement Field")
add_rect(gca())
savefig("../plots/rect_intro_lena_flow.png")

# rect disp
flow[100:130, 110:140]
showflow(flow[100:130, 110:140], figtitle="Displacement Field of Rectagle")
rect = patches.Rectangle((1,1),30,30,linewidth=1,edgecolor="r",facecolor="none")
ax = gca()
ax.add_patch(rect)
gcf()
savefig("../plots/rect_intro_lena_flow_zoom.png")


## lap const displacement
img, imgw, flow = gen_init(:lena, flow_args=[20, 300])

imgshow(img, figtitle="Target")
savefig("../plots/lap_const_lena_orig.png")

imgshow(imgw, figtitle="Source")
savefig("../plots/lap_const_lena_warped.png")

imgoverlay(img, imgw, figtitle="Blending of Target and Source")
savefig("../plots/lap_const_lena_bl_or_movement.png")

showflow(flow, figtitle="Constant Displacement field")
savefig("../plots/lap_const_lena_flow.png")

## Make images for testing

img, imgw, flow = gen_init(:lena, flow_args=[20, 60])

showflow(flow)

imgshow(imgw)

using Images, FileIO, ImageIO


save("../testimages/smoothly_varying/target.png", img)
save("../testimages/smoothly_varying/source.png", imgw)


## make tests images:
using Printf
#params:
num_tests = 3
basedir = "../testimages/"
folder = "settings_1/"
flow_args_smoothly_var = [20, 60]
flow_args_const = [20, 300]
img_select = :chess

if img_select == :chess
    folder_name = "chess_"
elseif img_select == :lena
    folder_name = "lena_"
end



for test in 1:num_tests
    name = folder_name * string(test) * "/"
    base_path = basedir*folder*name

    for setting in ["smoothly_varying", "constant"]
        if setting == "smoothly_varying"
            img, imgw, flow = gen_init(img_select, flow_args=flow_args_smoothly_var)
        elseif setting == "constant"
            img, imgw, flow = gen_init(img_select, flow_args=flow_args_const)
        end
        path = base_path * setting * "/"
        mkpath(path);

        # save the max displacement into a file
        f = open(path * "max_disp.txt", "w")
        magnitudes = map(x -> LAP_julia.vec_len(x), flow)
        max_mag = maximum(filter(!isnan, magnitudes))
        max_disp_str = @sprintf("%0.2f", max_mag)
        write(f, max_disp_str)
        close(f)

        save(path * "/target.png", img)
        save(path * "/source.png", imgw)
    end
end

## evaluate tests

methods = ["bunwarpj", "pflap_matlab"]

img, imgw, flow = gen_init(flow_args=[20, 120])
showflow(flow)

u_est, source_reg = pflap(img, imgw)

LAP_julia.mean(abs.(img .- imgw))
LAP_julia.mean(abs.(img .- source_reg))

dist_euclid = sqrt(sum((img .- imgw).^2))/length(img)
dist_euclid = sqrt(sum((img .- source_reg).^2))/length(img)

ssd(img, imgw)
ssd(img, source_reg)
ssd(img, source_reg_p)

imgshow(source_reg)

LAP_julia.mean(source_reg)


u_est_p, source_reg_p = sparse_pflap(img, imgw)

showflow(u_est_p)

imgshow(source_reg_p)
imgshow(img)

LAP_julia.mean(abs.(img .- source_reg_p))
dist_euclid = sqrt(sum(img .- imgw).^2)/length(img)


u_est, source_reg_p_ed, figs, Δ_u = sparse_pflap(img, imgw)

imgshow(source_reg_p_ed)
imgshow(img)


figs[1,5]


## single_lap


u_est = single_lap(img, imgw, 20, [41,41])

showflow(u_est)
LAP_julia.inpaint_nans!(u_est)
u_est_sm = LAP_julia.smooth_with_gaussian!(u_est, 20)

showflow(u_est_sm)


vec, edge = LAP_julia.gradient_points.gradient_magnitude(img)

imgshow(edge, figtitle="Edge Image")
savefig("../plots/sparse_lap_edge_image.png")

## sparse lap

img, imgw, flow = gen_init(:lena, flow_args=[10, 150])
showflow(flow, figtitle="Ground Truth Displacement")
savefig("../plots/sparse_lap_ground_truth_flow.png")

fhs = 15
window_size = [31, 31]

mask = parent(padarray(trues(size(img).-(2*fhs, 2*fhs)), Fill(false, (fhs, fhs), (fhs, fhs))))
inds = find_edge_points(img, mask=mask)
points = LAP_julia.inds_to_points(inds)
new_estim = single_lap_at_points(img, imgw, fhs, window_size, points, 3)
showflow(new_estim, disp_type=:sparse, figtitle="Estimated Displacement Vectors")
savefig("../plots/sparse_lap_estimate_image.png")

full_new_estim = interpolate_flow(new_estim, inds)
showflow(full_new_estim, figtitle="Interpolated Displacement Vectors")
savefig("../plots/sparse_lap_estimate_interpolated_image.png")


imgshow(edge, figtitle="Edge Image with Points")
PyPlot.scatter([ind[2] for ind in inds], [ind[1] for ind in inds], marker = :x); gcf()
savefig("../plots/sparse_lap_edge_image_with_points.png")
