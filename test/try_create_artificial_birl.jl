
using LAP_julia


view(img) = colorview(Gray, img);


img_size = (1234,1234)

img = gen_spaghetti(img_size)
view(img)

img, imgw, flow = gen_init(:spaghetti, img_args=[img_size], flow_args=[125])

view(img)
view(imgw)

timer = TimerOutput(string(args_dict["reg_alg"]))
method_kwargs = Dict(:display => false, :timer => timer, :match_source_histogram => true)
flow_est, source_reg, timer, time_in_secs = test_registration_alg(sparse_pflap_psnr,
                                                                  img,
                                                                  imgw,
                                                                  flow,
                                                                  timer=method_kwargs[:timer],
                                                                  method_kwargs=method_kwargs)

view(source_reg)
