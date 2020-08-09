
using LAP_julia: load_H, make_flow_from_H, gen_rand_points

## create a birl like structure for with the oxford affine dataset
base_path = "/Users/MrTrololord/Downloads/ox_affine/"
images_path = joinpath(base_path, "images")
folders = filter(x -> !isfile(x), readdir(base_path))



img = load(joinpath(base_path, folders[1], "img1.ppm"))


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



COUNT = 500

# create landmarks for target image --> img1
target_lnd_paths = []
for path in target_img_paths
    folder = split(path, "/")[end-1]
    target_img = load(path)
    target_lnd_path = joinpath(base_path, "landmarks", folder, "target.csv")
    target_landmarks = create_landmarks(target_img, COUNT, target_lnd_path)
    push!(target_lnd_paths, target_lnd_path)
    # target_lnd = create_landmarks()
    # lnd_path = joinpath(base_path, folder, )
end


function create_landmarks(img, count, save_path)
    inds = gen_rand_points(img, count, "Gridded")
    @info length(inds)
    df = DataFrame(Column1 = 1:length(inds), X = map(ind -> ind[2], inds), Y = map(ind -> ind[1], inds))
    CSV.write(save_path, df)
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
    LAP_julia.save_shift_landmarks(target_lnd_path, cur_flow, save_path)
end


# create loc table
tmp_df = CSV.read("/Volumes/davidkunz/Documents/bakalarka/birl_small/dataset_small_repaired.csv") |> DataFrame
names(df)
df = deepcopy(tmp_df)
delete!(df, 1:size(tmp_df,1))


i = 0
for folder in folders
    for k in 2:6
        global i = i + 1
        row = Dict()
        row[string("Column1")] = i
        # images
        row[string("Target image")] = joinpath(base_path, "images", folder, "img1.png")
        row[string("Source image")] = joinpath(base_path, "images", folder, "img$k.png")
        # lnds
        row[string("Target landmarks")] = joinpath(base_path, "landmarks", folder, "target.csv")
        row[string("Source landmarks")] = joinpath(base_path, "landmarks", folder, "source$k.csv")
        # sizes
        source_img = load(row["Source image"])
        row[string("Image size [pixels]")] = size(source_img)
        row[string("Image diagonal [pixels]")] = round(sqrt(sum(size(source_img).^2)), digits=1)
        row[string("status")] = "training"
        # empty values
        row[string("Warped target landmarks")] = Missing
        row[string("Warped source landmarks")] = Missing
        row[string("Execution time [minutes]")] = Missing

        push!(df, row)
    end
end

target_source_transform_paths


img = load(joinpath(base_path_bikes, "img1.png"))
k = 4
imgw = load(joinpath(base_path_bikes, "img$(k).png"))


brain_path = "/Users/MrTrololord/Google_Drive/cvut/bakalarka/cancer/AAPM-RT-MAC/RTMAC-LIVE-010/1.3.6.1.4.1.14519.5.2.1.1706.6003.102033049605502186962807064832/1.3.6.1.4.1.14519.5.2.1.1706.6003.218679427253271214825987989417/1-010.png"

brain = load(brain_path)

flow = gen_quad_flow(size(brain))
brainw = warp_img(brain, -real.(flow), -imag.(flow))

flow_est, source_reg = time_reg_alg(sparse_pflap, brain, brainw)

showflow(flow_est)
showflow(flow.*(-1), figtitle="truth")
