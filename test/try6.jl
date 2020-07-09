


function gen_flow(img_size)

    a, b, c = randn(), randn(), randn()
    f(z) = a + b*z + c*z^2 # + d*z^3 + e*z^4

    X = ones(img_size[1]) * collect(range(0,1,length=img_size[2]))'
    Y = collect(range(0,1,length=img_size[1])) * ones(img_size[2])'
    A = Y .* im + X

    B = f.(A)
    max_len = maximum(LAP_julia.vec_len.(B))
    B = B./max_len
    return B
end

B = gen_flow((12,12)); showflow(B)

println(B)

f(1 * im)

imgshow(rand(12,123))
