#-------------------------------------------------------------------------------
# Set of useful functions related to LinearModel struct.
#-------------------------------------------------------------------------------

struct LinearModel

    imin  :: Int64
    c     :: Float64
    g     :: AbstractVector
    xbase :: AbstractVector
    fval  :: AbstractVector
    dst   :: AbstractVector
    Y     :: AbstractMatrix

end

LinearModel(imin, xbase, fval, Y) = begin

    m, n = size(Y)
    dst = zeros(Float64, m)

    for j =1:n
        for i = 2:m
            Y[i, j] -= xbase[j]
            dst[i - 1] += Y[i, j] ^ 2.0
        end
    end
    @. dst = sqrt(dst)
    
    g = A \ (fval[2:end] .- fval[1])
    c = fval[1] - dot(g, xbase)

    LinearModel(imin, c, g, xbase, fval, dst, Y)

end

