
using LAP_julia: load_H, make_flow_from_H, gen_rand_points, get_valid_landmarks
using DataFrames, CSV

## create a birl like structure for with the oxford affine dataset
base_path = "/Users/MrTrololord/Downloads/head_mri_homography/"
target_images_path = joinpath(base_path, "images", "target")
target_images_fnames = filter(x -> !isfile(x), readdir(target_images_path))
target_images_paths = map(fname -> joinpath(target_images_path, fname), target_images_fnames)

COUNT = 500
max_disp = 40
# flow_gen_func(img_size) = gen_quad_flow(img_size, max_disp)
flow_gen_func(img_size) = gen_homo_flow(img_size, max_disp)

for target_img_path in target_images_paths
    k = parse(Int, target_img_path[findfirst(r"[0-9]+", target_img_path)])
    target_img = load(target_img_path)
    # create save paths
    target_land_path = joinpath(base_path, "landmarks", "target", "landmarks$(k).csv")
    source_land_path = joinpath(base_path, "landmarks", "source", "landmarks$(k).csv")
    source_img_path = joinpath(base_path, "images", "source", "img$(k).png")
    # generate rand inds in target
    target_inds = gen_rand_points(target_img, COUNT, "Gridded")
    cur_flow = flow_gen_func(size(target_img))
    # edit target inds and save
    target_inds = LAP_julia.get_valid_landmarks(cur_flow, target_inds)
    LAP_julia.save_landmarks(target_inds, target_land_path)
    # warp target inds and save
    source_points = LAP_julia.move_landmarks(transpose(LAP_julia.inds_to_points(target_inds)), cur_flow
    # rng = extrema.(indices_spatial(cur_flow))
    # @info source_points[:,1] size(cur_flow) rng
    # @assert all(LAP_julia.is_in_bounds.(source_points[:,1], rng[1]...) .& LAP_julia.is_in_bounds.(source_points[:,2], rng[2]...))
    LAP_julia.save_landmarks(source_points, source_land_path)
    # LAP_julia.save_shift_landmarks(target_land_path, cur_flow.*(-1), source_land_path)
    # create source img
    source = warp_img(target_img, -real.(cur_flow), -imag.(cur_flow))
    save(source_img_path, source)
    @info k target_img_path size(target_inds)
end


# create loc table
tmp_df = CSV.read("/Volumes/davidkunz/Documents/bakalarka/birl_small/dataset_small_repaired.csv") |> DataFrame

df = deepcopy(tmp_df)
delete!(df, 1:size(tmp_df,1))

# rename!(df, Symbol.(replace.(string.(names(df)), Ref(" "=> "_"))))
# rename!(df, Symbol.(names(df)))

names(df)
eltypes(df)

abc = [:a , :b, :c, :d, :e, :f, :g, :h, :i, :j, :k]
new_df = DataFrame(a = eltypes(df)[1][],
                   b = eltypes(df)[2][],
                   c = eltypes(df)[3][],
                   d = eltypes(df)[4][],
                   e = eltypes(df)[5][],
                   f = eltypes(df)[6][],
                   g = eltypes(df)[7][],
                   h = eltypes(df)[8][],
                   i = eltypes(df)[9][],
                   j = eltypes(df)[10][],
                   k = eltypes(df)[11][])

function rename_dict_keys(from_name, to_name, dict)
    new_dict = Dict()
    rename_dict = Dict(from_name .=> to_name)
    for (key, value) in rename_dict
        new_dict[value] = dict[key]
    end
    return new_dict
end

experiment_count = count(x -> !isfile(x), readdir(target_images_path))


for k in 1:experiment_count
    row = Dict{String}{Any}()
    row["Column1"] = k
    # images
    row["Target image"] = joinpath("images", "target", "img$k.png")
    row["Source image"] = joinpath("images", "source", "img$k.png")
    # lnds
    row["Target landmarks"] = joinpath("landmarks", "target", "landmarks$k.csv")
    row["Source landmarks"] = joinpath("landmarks", "source", "landmarks$k.csv")
    # sizes
    source_img = load(joinpath(base_path, "images", "source", "img$k.png"))
    row["Image size [pixels]"] = string(size(source_img))
    row["Image diagonal [pixels]"] = round(sqrt(sum(size(source_img).^2)), digits=1)
    row["status"] = "training"
    # empty values
    row["Warped target landmarks"] = missing
    row["Warped source landmarks"] = missing
    row["Execution time [minutes]"] = missing

    @info rename_dict_keys(names(tmp_df), abc, row)
    push!(new_df, rename_dict_keys(names(tmp_df), abc, row))
end

new_df

# rename back
renamed = rename!(new_df, names(tmp_df))

CSV.write(joinpath(base_path, "loc_table.csv"), renamed)
