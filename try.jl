using LAP_julia, TimerOutputs

for _ in 1:3
    img, imgw, flow = gen_init(:lena, :tiled, flow_args=[15, 30])

    println("ONCE")
    timer=TimerOutput("pf lap");
    method_kwargs =Dict(:timer => timer, :display => false, :max_repeats => 1)
    timer = @time flow_est, source_reg, timer, results = test_registration_alg(pflap, img, imgw, flow, [], method_kwargs, timer=timer, display=true);
end
