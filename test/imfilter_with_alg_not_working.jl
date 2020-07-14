using ImageFiltering


img = rand(50,50)
alg = Algorithm.FIR
factor = ImageFiltering.centered(ones(13))
factored_ker = kernelfactors((factor, factor))

# works
imfilter(img, factored_ker, "symmetric")

# breaks
imfilter(img, factored_ker, "symmetric", alg)
