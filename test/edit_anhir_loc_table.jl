
using CSV, FileIO, DataFrames

base_path = "/Users/MrTrololord/Google_Drive/cvut/bakalarka/anhir/"
base_path = "/Volumes/davidkunz/Documents/bakalarka/birl_small/"

loc_table = CSV.read(base_path * "location_table.csv") |> DataFrame
loc_table = CSV.read(base_path * "dataset_small.csv") |> DataFrame

train_rows = loc_table[loc_table[:status] .== "training", :]

non_mammary_gland2 = train_rows[map(x-> !occursin("mammary-gland_2", x), collect(train_rows["Source image"])), :]


train_rows = non_mammary_gland2
# train_rows[map(x-> !occursin("mammary-gland_2", x), collect(train_rows["Source image"])), :]
# map(x-> occursin("mammary-gland_2", x), collect(train_rows["Source image"]))

#first k rows:
k = 15
subset = train_rows[1:k, :]

image_folder = "images/"
landmark_folder = "landmarks/"

for row in eachrow(train_rows)
    row["Source image"] = image_folder * row["Source image"]
    row["Target image"] = image_folder * row["Target image"]
    row["Source landmarks"] = landmark_folder * row["Source landmarks"]
    row["Target landmarks"] = landmark_folder * row["Target landmarks"]
end

show(train_rows)



CSV.write(joinpath(base_path, "dataset_small_repaired.csv"), train_rows)


##
# train_rows[1, "Source landmarks"]
landmark_path = joinpath(base_path, train_rows[1, "Source landmarks"])

source_landmarks = CSV.read(landmark_path) |> DataFrame

locations = [source_landmarks[:, "Y"] source_landmarks[:, "X"]]

rand_flow = gen_quad_flow((2762, 2099), 12)

LAP_julia.move_landmarks(locations, rand_flow)

itp = Interpolations.interpolate(rand_flow, BSpline(Linear()))

loc = locations[1,:]
itp(loc...)

shifted_locations = similar(locations)
[real(itp(locations[1,:]...)), imag(itp(locations[1,:]...))]

for k in 1:size(locations, 1)
    shifted_locations[k,:] .= locations[k,:] .+ [real(itp(locations[k,:]...)), imag(itp(locations[k,:]...))]
end

shifted_locations

showflow(rand_flow)

addpoints(locations)
addpoints(shifted_locations)

img, imgw, flow = gen_init();
locs = LAP_julia.inds_to_points(rand(CartesianIndices(size(flow)), 60))'

shifted_locs = LAP_julia.move_landmarks(locs, flow)

loc = CartesianIndex(2,123)
addpoints([loc])

showflow(flow)
addpoints(locs)
addpoints(shifted_locs)


## save locs

locations = [source_landmarks[:, "Y"] source_landmarks[:, "X"]]
shifted_landmark_df = DataFrame(Column1 = 0:(size(locations,1)-1), X = locations[:,2], Y = locations[:,1])
shifted_landmark_df == source_landmarks

save_path = joinpath(args_dict["output_path"], args_dict["land_warped_fname"])
CSV.write(save_path, shifted_landmark_df)


function shift_landmarks(landmarks_path, flow, args_dict)
    source_landmarks_df = CSV.read(landmark_path) |> DataFrame
    locations = [source_landmarks_df[:, "Y"] source_landmarks_df[:, "X"]]
    shifted_landmarks = LAP_julia.move_landmarks(locations, flow)
    shifted_landmarks_df = DataFrame(Column1 = 0:(size(locations,1)-1), X = shifted_landmarks[:,2], Y = shifted_landmarks[:,1])
    save_path = joinpath(args_dict["output_path"], args_dict["land_warped_fname"])
    CSV.write(save_path, shifted_landmark_df)
end
