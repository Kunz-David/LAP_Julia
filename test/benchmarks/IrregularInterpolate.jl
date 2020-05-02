module IrregularInterpolate

using LinearAlgebra

export PolyharmonicSpline, interpolate


struct PolyharmonicSpline{T}
    dim::Int
    order::Int
    coeff::Vector{T}
    centers::Matrix{T}
    error::T
end

function polyharmonicK(r, K)
    if iseven(K)
        return r < 1 ? r^(K-1) * log(r^r) : (r^K) * log(r)
    else
        return r^K
    end
end


function PolyharmonicSpline(K::Int, centers::Matrix{T}, values::Array{T}; s = zero(T)) where {T}
    N, dim = size(centers)
    N != length(values) && throw(DimensionMismatch())

    phK = Matrix{T}(undef, N, N)
    for i in 1:N
        for j in 1:N
            dist = 0.0
            for d in 1:dim
                dist += (centers[i, d] - centers[j, d])^2
            end
            phK[i, j] = polyharmonicK(sqrt(dist), K)
        end
    end

    A = copy(phK)
    B = zeros(T, N, dim+1)
    for i in 1:N
        B[i, 1] = 1
        B[i, 2:end] = centers[i, :]
    end

    A .+= s .* Diagonal(I, N)

    L = [A B; B' zeros(T, dim+1, dim+1)]
    w = pinv(L) * [values; zeros(T, dim+1)]

    ivalues = zeros(T, N)
    for i in 1:N
        tmp = 0.0
        for j in 1:N
            tmp += w[j] * phK[i, j]
        end
        tmp += w[N+1]
        for j in 2:dim+1
            tmp += w[N+j] * centers[i, j-1]
        end
        ivalues[i] = tmp
    end
    error = norm(values .- ivalues)

    return PolyharmonicSpline(dim, K, w, centers, error)
end

function PolyharmonicSpline(K::Int, centers::Vector{T}, values::Vector{T}; s = zero(T)) where {T}
    # PolyharmonicSpline(K, centers'', values, s = s)
    PolyharmonicSpline(K, centers, values, s = s)
end

function interpolate(S::PolyharmonicSpline{T}, x::Matrix{T}) where {T}
    m, n = size(x)
    n != S.dim && throw(DimensionMismatch("$m != $(S.dim)"))

    l = length(S.coeff) - (n+1)

    interpolates = zeros(T, m)
    for i in 1:m
        tmp = 0.0
        for j in 1:l
            dist = 0.0
            for d in 1:n
                dist += (x[i, d] - S.centers[j, d])^2
            end
            tmp += S.coeff[j] * polyharmonicK(sqrt(dist), S.order)
        end
        tmp += S.coeff[l+1]
        for j in 2:n+1
            tmp += S.coeff[l+j] * x[i, j-1]
        end
        interpolates[i] = tmp
    end
    return interpolates
end

function interpolate(S::PolyharmonicSpline{T}, x::Vector{T}) where {T}
    return interpolate(S, x)
end

function interpolate(S::PolyharmonicSpline{T}, x::Vector{T}, y::Vector{T}) where {T}
    return interpolate(S, [x y])
end

function interpolate(S::PolyharmonicSpline{T}, x::Vector{T}, y::Vector{T}, z::Vector{T}) where {T}
    return interpolate(S, [x y z])
end

end  # module
