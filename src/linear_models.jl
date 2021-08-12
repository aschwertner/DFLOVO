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

function update!(model::LinearModel, t, fval_d, d)

    # Updates Y set.
    for i = 1:length(model.dst)
        @. model.Y[i, :] -= d
    end

    # Puts the old xbase in the t-th vector in Y.
    @. model.Y[t, :] =  - d
    model.fval[t + 1] = model.fval[1]

    # Updates the distance vector.
    for i = 1:len(model.dst)
        model.dst[i] = norm(model.Y[i, :])
    end

    # Adds the new center
    model.fval[1] = fval_d
    @. model.xbase += d 

    return rebuild!(model)

end