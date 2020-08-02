
using LAP_julia


view(img) = colorview(Gray, img);


img_size = (124,121)

img, imgw, flow = gen_init(:spaghetti, img_args=[img_size], flow_args=[25])

view(img)
view(imgw)

timer = TimerOutput(string(args_dict["reg_alg"]))
method_kwargs = Dict(:display => true, :timer => timer, :match_source_histogram => true)
flow_est, source_reg, timer, time_in_secs = test_registration_alg(sparse_pflap_psnr,
                                                                  img,
                                                                  imgw,
                                                                  flow,
                                                                  timer=method_kwargs[:timer],
                                                                  method_kwargs=method_kwargs)


flow_est
LAP_julia.vec_len.(flow_est)
showflow(flow_est)
view(source_reg)
view(img)
