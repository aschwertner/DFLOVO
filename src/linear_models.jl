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

function LinearModel(imin, xbase, fval, Y)

    m, n = size(Y)
    dst = zeros(Float64, m)

    for j =1:n
        for i = 2:m
            Y[i, j] -= xbase[j]
            dst[i - 1] += Y[i, j] ^ 2.0
        end
    end
    @. dst = sqrt(dst)
    
    g = A \ ( fval[2:end] .- fval[1] )
    c = fval[1] - dot( g, xbase )

    return LinearModel(imin, c, g, xbase, fval, dst, Y)

end

( model::LinearModel )( x::AbstractVector ) = model.c + dot( model.g, x )

∇( model::LinearModel )( x::AbstractVector ) = model.g

∇2( model::LinearModel )( x::AbstractVector ) = false * I

