
# julia test/try7.jl /Users/MrTrololord/Google_Drive/cvut/bakalarka/BIRL/data-images/images/artificial_reference.jpg /Users/MrTrololord/Google_Drive/cvut/bakalarka/BIRL/data-images/images/artificial_moving-affine.jpg /Users/MrTrololord/Google_Drive/cvut/bakalarka/BIRL/data-images/landmarks/artificial_moving-affine.csv  /Users/MrTrololord/Google_Drive/cvut/bakalarka/LAP_julia/test/outputs/ sparse_pflap

function move_landmarks(locations, flow)
    return map(loc -> loc .+ (real.(flow)[loc...], imag.(flow)[loc...]), locations)
end


flow[1,1]
move_landmarks([(1,1)], flow)

args_dict = Dict()



target = load_image_gray("/Users/MrTrololord/Google_Drive/cvut/bakalarka/BIRL/data-images/images/artificial_reference.jpg")
source = load_image_gray("/Users/MrTrololord/Google_Drive/cvut/bakalarka/BIRL/data-images/images/artificial_moving-affine.jpg")

args_dict["reg_alg"] = sparse_pflap_psnr
timer = TimerOutput(string(args_dict["reg_alg"]));

showflow(flow)

# target, source, flow = gen_init(:spaghetti, img_args=[(12,10)], flow_args=[2])


method_kwargs =Dict(:timer => timer, :display => true, :max_repeats => 3, :match_source_histogram => false)
flow_est, source_reg, timer, time_in_secs, (figs,) = time_reg_alg(args_dict["reg_alg"],
                                                         target,
                                                         source,
                                                         timer=timer,
                                                         method_kwargs=method_kwargs);
println()

showflow(flow_est)
imgoverlay(target, source_reg)
imgoverlay(target, source)


figs[1,1,3]


target, source, flow = gen_init()

args_dict["reg_alg"] = pflap
timer = TimerOutput(string(args_dict["reg_alg"]))
method_kwargs = Dict(:timer => timer, :display => false, :max_repeats => 2, :match_source_histogram => false)
flow_est, source_reg, timer, time_in_secs = time_reg_alg(args_dict["reg_alg"],
                                                         target,
                                                         source,
                                                         timer=timer,
                                                         method_kwargs=Dict(:display => true, :timer => timer))
println()

imgshow(source_reg)
imgshow(target)
imgoverlay(target, source_reg)


#
# y = [1 2; 3 4]
# kern = centered([0 1 0;1 -4 1;0 1 0])
#
# imfilter(y, kern, Fill(0,kern))
#
# import LAP_julia: estimation_noise_variance
#
# LAP_julia.estimation_noise_variance(target)
#
# img, imgw, flow = gen_init()
#
# showflow(flow)
# imgshow(imgw)
#
# LAP_julia.mapped_out(flow.*4)
#
# imgw_copy = deepcopy(imgw)
# imgw_copy[LAP_julia.mapped_out(flow)] .= 0
# imgshow(imgw_copy)
# imgshow(imgw, origin_left_bot = true)
#


assess_psnr(target, source)

assess_psnr(img, imgw)
assess_psnr(img, source_reg)

imgshow(source_reg)

imgshow(img)
imgshow(imgw)
imgoverlay(img,imgw)
imgoverlay(source_reg, img)

timer=TimerOutput("pf lap");
method_kwargs =Dict(:timer => timer, :display => false, :max_repeats => 3, :match_source_histogram => false)
flow_est, source_reg, timer, results = test_registration_alg(pflap, img, imgw, flow, method_kwargs=method_kwargs, timer=timer);


showflow(flow, figtitle="truth flow")

figs[1,1,1]
figs[6,1,1]


img, imgw, flow = gen_init(flow_args=[30]);

timer=TimerOutput("pflap reg");
method_kwargs =Dict(:timer => timer, :display => false, :max_repeats => 3, :match_source_histogram => false)
flow_est, source_reg, timer, results = test_registration_alg(pflap, img, imgw, flow, method_kwargs=method_kwargs, timer=timer);


timer=TimerOutput("sparse pflap reg"); # fix
method_kwargs =Dict(:timer => timer, :display => true, :max_repeats => 3, :match_source_histogram => false)
flow_est, source_reg, timer, results = test_registration_alg(sparse_pflap, img, imgw, flow, method_kwargs=method_kwargs, timer=timer);

timer=TimerOutput("sparse pflap psnr reg"); # fix
method_kwargs =Dict(:timer => timer, :display => true, :max_repeats => 3, :match_source_histogram => false)
flow_est, source_reg, timer, results, (figs,) = test_registration_alg(sparse_pflap_psnr, img, imgw, flow, method_kwargs=method_kwargs, timer=timer);

timer=TimerOutput("sparse lap reg");
method_kwargs =Dict(:timer => timer, :display => false)
flow_est, source_reg, timer, results = test_registration_alg(sparse_lap, img, imgw, flow, method_args=[30], method_kwargs=method_kwargs, timer=timer);

timer=TimerOutput("lap reg");
method_kwargs =Dict(:timer => timer, :display => false)
flow_est, source_reg, timer, results = test_registration_alg(lap, img, imgw, flow, method_args=[30], method_kwargs=method_kwargs, timer=timer);



imgshow(img, origin_left_bot = true, figtitle="target")
imgshow(imgw, origin_left_bot = true, figtitle="source")
imgshow(source_reg, origin_left_bot = true, figtitle="source_reg")

imgshow(img .- source_reg, origin_left_bot = true, figtitle="diff")

showflow(flow.*(-1), figtitle="target flow")
showflow(flow_est, figtitle="estimated flow")
showflow((flow_est .+ flow)[30:end-30, 30:end-30], figtitle="flow diff")


###

#anhir loc table subset

base_path = "/Users/MrTrololord/Google_Drive/cvut/bakalarka/anhir/"
loc_table = CSV.read(base_path * "location_table.csv")

train_rows = loc_table[loc_table[:status] .== "training", :]

#first k rows:
k = 15
subset = train_rows[1:k, :]

CSV.write(joinpath(base_path, "subset_table.csv"), subset)



##
using LAP_julia, TimerOutputs

target, source = gen_anhir(mutate=true, diag_pixels=500)

imgshow(target)
imgshow(source)

timer=TimerOutput("pflap reg");
method_kwargs =Dict(:timer => timer, :display => false)
flow_est, source_reg, timer, results = time_reg_alg(pflap, target, source, method_kwargs=method_kwargs, timer=timer);


timer=TimerOutput("sparse pflap reg");
method_kwargs =Dict(:timer => timer, :display => false)
flow_est, source_reg, timer, results = time_reg_alg(sparse_pflap_psnr, target, source, method_kwargs=method_kwargs, timer=timer);


imgshow(source_reg)

showflow(flow_est)

A = rand(1234,1234)

@benchmark LAP_julia.resize_to_diag_size(A, 500)


flow_est























## tmp

args_dict = Dict("img_path" => "/Users/MrTrololord/Google_Drive/cvut/bakalarka/BIRL/data-images/images/artificial_reference.jpg",
                 "imgw_path" => "/Users/MrTrololord/Google_Drive/cvut/bakalarka/BIRL/data-images/images/artificial_moving-affine.jpg",
                 "moving_landmark_path" => "/Users/MrTrololord/Google_Drive/cvut/bakalarka/BIRL/data-images/landmarks/artificial_moving-affine.csv",
                 "output_path" => "/Users/MrTrololord/Google_Drive/cvut/bakalarka/LAP_julia/test/outputs/",
                 "reg_alg" => "sparse_pflap_psnr",
                 "time_fname" => "exec_time_in_secs.txt",
                 "land_warped_fname" => "warped_landmarks.csv",
                 "diag_pix" => 500)


base_path = "mount/birl_small/"
args_dict = Dict("img_path" => joinpath(base_path, "images/COAD_01/scale-5pc/S3.jpg"),
                 "imgw_path" => joinpath(base_path, "images/COAD_01/scale-5pc/S1.jpg"),
                 "moving_landmark_path" => joinpath(base_path, "landmarks/COAD_01/scale-5pc/S3.csv"),
                 "output_path" => joinpath(base_path, "tmp/"),
                 "reg_alg" => "sparse_pflap_psnr",
                 "time_fname" => "exec_time_in_secs.txt",
                 "land_warped_fname" => "warped_landmarks.csv",
                 "diag_pix" => 500)




println("Running with the arguments:")

println("******ARGS:********")
args_dict["reg_alg"] = getfield(LAP_julia, Symbol(args_dict["reg_alg"]))

target = load_image_gray(args_dict["img_path"])
source = load_image_gray(args_dict["imgw_path"])

resized_target, target_resize_ratio = LAP_julia.resize_to_diag_size(target, args_dict["diag_pix"])
resized_source, source_resize_ratio = LAP_julia.resize_to_diag_size(source, args_dict["diag_pix"])

ready_target, ready_source = LAP_julia.pad_images(resized_target, resized_source)

# Take the 3 run of the reg alg, allowing for julia compilation.
for _ in 1:3
    timer = TimerOutput(string(args_dict["reg_alg"]))
    flow_est, source_reg, timer, time_in_secs = time_reg_alg(args_dict["reg_alg"],
                                                             ready_target,
                                                             ready_source,
                                                             timer=timer,
                                                             method_kwargs=Dict(:display => false, :timer => timer))
    global flow_est, source_reg, timer, time_in_secs
    println()
end

# save time
save_time(time_in_secs, args_dict)

# rescale flow_est
flow_est = imresize(flow_est.*(1/source_resize_ratio), ratio=1/source_resize_ratio)

# move landmarks
save_shift_landmarks(args_dict["moving_landmark_path"], flow_est, args_dict)

source_landmarks_df = CSV.read(args_dict["moving_landmark_path"]) |> DataFrame
locations = [source_landmarks_df[:, "Y"] source_landmarks_df[:, "X"]]
shifted_landmarks = LAP_julia.move_landmarks(locations, flow_est)

addpoints(locations)
addpoints(shifted_landmarks)
