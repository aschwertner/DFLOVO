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

    n = length(xbase)
    dst = zeros(Float64, n - 1)

    #for j =1:n
    #    for i = 1:n
    #        Y[i, j] -= xbase[j]
    #        dst[i - 1] += Y[i, j] ^ 2.0
    #    end
    #end
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

function fromscratch!(model::LinearModel, func_list, δ, a, b)

    n = length(model.dst)

    # Sets Y to a empty matrix
    @. model.Y = 0.0

    for i = 1:n
        
        # Computes the new interpolation points and the distances
        α = min( u[i] - model.xbase[i], δ )
        model.Y[i, i] = α
        model.dst[i] = α

        # Computes the function values
        model.xbase[i] += α
        model.fval[i + 1] = func_list[imin](model.xbase)
        model.xbase[i] -= α

    end

    return rebuild!(model)

end

function rebuild_model!(model::LinearModel)

    # Computes the new c and g
    # !!! Very inefficient !!!

    model.g = model.Y \ ( model.fval[2:end] .- model.fval[1] )
    model.c = model.fval[1] - dot( model.g, model.xbase )

    return LinearModel(model.imin, model.c, model.g, model.xbase, model.fval, model.dst, model.Y )

end