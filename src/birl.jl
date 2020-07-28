using FileIO, Images, Colors

function register(alg, target_path, source_path)

    target, source = load_anhir_image_pair(target_path, source_path, base_path="")



    method_args = []
    method_kwargs = Dict()

    flow_est, source_reg = alg(target, source)

end

img, imgw, flow = gen_init(:chess, img_args=[50, 10])

point_count = 700
spacing = 100

full_flow_estim, source_reg, flow_estim_at_inds, inds = sparse_lap(img, imgw, 25, point_count=point_count, spacing=spacing)
showflow(create_sparse_flow_from_full(full_flow_estim, inds)); addpoints(inds);
println(length(inds))
