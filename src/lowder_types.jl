#-------------------------------------------------------------------------------
# LOWDER types
#-------------------------------------------------------------------------------

abstract type AbstractModel end

struct LinearModel <: AbstractModel

    n     :: Int64              # Model dimension.
    m     :: Int64              # Number of interpolation points.
    imin  :: Ref{Int64}         # Index of the selected function in 'func_list'.
    kopt  :: Ref{Int64}         # Index of the best point so far in 'Y'.
    c     :: Ref{Float64}       # Constant of the model.
    g     :: AbstractVector     # Gradient of the model.
    xbase :: AbstractVector     # Origin of the sample set.
    xopt  :: AbstractVector     # Best point so far (in terms of function values).
    fval  :: AbstractVector     # Set of the function values of the interpolation points.
    dst   :: AbstractVector     # Distances between 'xbase' and other interpolation points.
    Y     :: AbstractMatrix     # Set of interpolation points, shifted from the center of the sample set 'xbase'.

end