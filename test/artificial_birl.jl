
using LAP_julia: load_H, make_flow_from_H, gen_rand_points
using DataFrames, CSV

## create a birl like structure for with the oxford affine dataset
base_path = "/Users/MrTrololord/Downloads/head_mri/"
target_images_path = joinpath(base_path, "images", "target")
target_images_fnames = filter(x -> !isfile(x), readdir(target_images_path))
target_images_paths = map(fname -> joinpath(target_images_path, fname), target_images_fnames)

COUNT = 500
max_disp = 50
flow_gen_func(img_size) = gen_quad_flow(img_size, max_disp)

for (k, target_img_path) in enumerate(target_images_paths)
    target_img = load(target_img_path)
    # create save paths
    target_land_path = joinpath(base_path, "landmarks", "target", "landmarks$(k).csv")
    source_land_path = joinpath(base_path, "landmarks", "source", "landmarks$(k).csv")
    source_img_path = joinpath(base_path, "images", "source", "img$(k).png")
    # generate rand inds in target
    target_inds = gen_rand_points(target_img, COUNT, "Gridded")
    cur_flow = flow_gen_func(size(img))
    # edit target inds and save
    target_inds = get_valid_landmarks(cur_flow, target_inds)
    LAP_julia.save_landmarks(target_inds, target_land_path)
    # warp target inds and save
    LAP_julia.save_shift_landmarks(target_land_path, cur_flow.*(-1), source_land_path)
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
