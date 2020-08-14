

base_path = "/Volumes/davidkunz/Documents/bakalarka/mri_head_homography"
target_images_path = joinpath(base_path, "images", "target")
target_images_fnames = filter(x -> !isfile(x), readdir(target_images_path))
target_images_paths = map(fname -> joinpath(target_images_path, fname), target_images_fnames)


for target_img_path in target_images_paths
    k = parse(Int, target_img_path[findfirst(r"[0-9]+", target_img_path)])
    #path
    source_img_path = joinpath(base_path, "images", "source", "img$(k).png")
    target = load(target_img_path)
    source = load(source_img_path)

    source_jpg_path = source_img_path[1:end-3]*"jpg"
    target_jpg_path = target_img_path[1:end-3]*"jpg"
    @info k source_jpg_path target_jpg_path
    # resave
    save(target_jpg_path, target)
    save(source_jpg_path, source)
end

source_images_path = joinpath(base_path, "images", "source")

for img in filter(x -> occursin(".png", x), readdir(source_images_path))
    run(`rm $(joinpath(source_images_path, img))`)
end

## remane to jpg

tmp_df = CSV.read("/Volumes/davidkunz/Documents/bakalarka/mri_head_homography/loc_table_png.csv") |> DataFrame

names(tmp_df)

collect(tmp_df["Source image"])
new_names_source = map(x -> replace(x, "png" => "jpg", count=1), collect(tmp_df["Source image"]))
new_names_target = map(x -> replace(x, "png" => "jpg", count=1), collect(tmp_df["Target image"]))

df = deepcopy(tmp_df)
df["Source image"] = new_names_source
df["Target image"] = new_names_target

CSV.write("/Volumes/davidkunz/Documents/bakalarka/mri_head_homography/loc_table.csv", df)
