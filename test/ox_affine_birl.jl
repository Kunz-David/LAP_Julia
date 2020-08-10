
using LAP_julia: load_H, make_flow_from_H, gen_rand_points

## create a birl like structure for with the oxford affine dataset
base_path = "/Users/MrTrololord/Downloads/ox_affine/"
images_path = joinpath(base_path, "images")
folders = filter(x -> !isfile(x), readdir(images_path))



#converted to png
for folder in folders
    my_path = joinpath(base_path, folder)
    for file in filter(x -> occursin("pgm", x), readdir(my_path))
        @info file
        img_path = joinpath(my_path, file)
        file_name = split(file, ".")[1]
        out_path = joinpath(my_path, file_name*".png")
        run(`convert $img_path $out_path`)
    end
end


# all pngs
image_paths = []
for folder in folders
    my_path = joinpath(images_path, folder)
    @info my_path
    for file in filter(x -> occursin("png", x), readdir(my_path))
        push!(image_paths, joinpath(my_path, file))
    end
end

# delete png base path
deleting = "/Users/MrTrololord/Downloads/"
image_paths_relative = map(x -> replace(x, deleting => ""), image_paths)


# get target images --> img1
target_img_paths = []
for folder in folders
    my_path = joinpath(images_path, folder)
    for file in filter(x -> occursin("img1.png", x), readdir(my_path))
        push!(target_img_paths, joinpath(my_path, file))
    end
end

target_img_paths


# paths of target_source_transform triads
target_source_transform_paths = []
for folder in folders
    my_path = joinpath(images_path, folder)
    for target in filter(x -> occursin("img1.png", x), readdir(my_path))
        target_path = joinpath(my_path, target)
        for source in filter(x -> !occursin("img1.png", x) && occursin("png", x), readdir(my_path))
            source_path = joinpath(my_path, source)
            img_number = parse(Int, source[end-4])
            tranform_file = filter(x -> !occursin("img", x) && occursin(string(img_number), x), readdir(my_path))[1]
            tranform_path = joinpath(my_path, tranform_file)
            push!(target_source_transform_paths, [target_path, source_path, tranform_path])
        end
    end
end
target_source_transform_paths = permutedims(hcat(target_source_transform_paths...))

target_source_transform_paths_relative = map(x -> replace(x, deleting => ""), target_source_transform_paths)


using DataFrames, CSV
COUNT = 500


# filter target landmarks not it source img
# create target_landmarks
for folder in folders
    target_img = load(joinpath(base_path, "images", folder, "img1.png"))
    target_lnd_path = joinpath(base_path, "landmarks", folder, "target.csv")
    target_inds = gen_rand_points(target_img, COUNT, "Gridded")
    for k in 2:6
        # get flow
        transform_path = joinpath(base_path, "images", folder, "H1to$(k)p")
        cur_H = LAP_julia.load_H(transform_path)
        cur_flow = LAP_julia.make_flow_from_H(cur_H, size(target_img))
        # shift and filter landmarks
        target_inds = LAP_julia.get_valid_landmarks(cur_flow.*(-1), target_inds)
    end
    @info size(target_inds)
    @info size(target_img)
    LAP_julia.save_landmarks(target_inds, target_lnd_path)
end


# warp target landmarks to source. and create
for (k, (target_path, source_path, transform_path)) in enumerate(eachrow(target_source_transform_paths))
    @info (k, target_path, source_path, transform_path)
    folder = split(target_path, "/")[7]
    target_lnd_path = joinpath(base_path, "landmarks", folder, "target.csv")
    cur_H = LAP_julia.load_H(transform_path)
    target_img = load(target_path)
    cur_flow = LAP_julia.make_flow_from_H(cur_H, size(target_img))
    @info size(cur_flow)
    img_number = transform_path[end-1]
    save_path = joinpath(base_path, "landmarks", folder, "source$img_number.csv")
    @info save_path
    LAP_julia.save_shift_landmarks(target_lnd_path, cur_flow.*(-1), save_path)
end


# create loc table
tmp_df = CSV.read("/Volumes/davidkunz/Documents/bakalarka/birl_small/dataset_small_repaired.csv") |> DataFrame

df = deepcopy(tmp_df)
delete!(df, 1:size(tmp_df,1))

rename!(df, Symbol.(replace.(string.(names(df)), Ref(" "=> "_"))))
rename!(df, Symbol.(names(df)))

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

i = 0
for folder in folders
    for k in 2:6
        global i = i + 1
        row = Dict{String}{Any}()
        row["Column1"] = i
        # images
        row["Target image"] = joinpath("images", folder, "img1.png")
        row["Source image"] = joinpath("images", folder, "img$k.png")
        # lnds
        row["Target landmarks"] = joinpath("landmarks", folder, "target.csv")
        row["Source landmarks"] = joinpath("landmarks", folder, "source$k.csv")
        # sizes
        source_img = load(joinpath(base_path, "images", folder, "img$k.png"))
        row["Image size [pixels]"] = string(size(source_img))
        row["Image diagonal [pixels]"] = round(sqrt(sum(size(source_img).^2)), digits=1)
        row["status"] = "training"
        # empty values
        row["Warped target landmarks"] = missing
        row["Warped source landmarks"] = missing
        row["Execution time [minutes]"] = missing
        global row

        @info rename_dict_keys(names(tmp_df), abc, row)
        push!(new_df, rename_dict_keys(names(tmp_df), abc, row))
    end
end

new_df

# rename back
renamed = rename!(new_df, names(tmp_df))

CSV.write(joinpath(base_path, "loc_table.csv"), renamed)
