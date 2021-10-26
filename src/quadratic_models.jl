#-------------------------------------------------------------------------------
# Set of useful functions related to LinearModel struct.
#-------------------------------------------------------------------------------

export QuadraticModel

struct QuadraticModel <: AbstractModel

    n     :: Int64              # Model dimension.
    m     :: Int64              # Number of interpolation points.
    imin  :: Ref{Int64}         # Index of the selected function in 'func_list'.
    kopt  :: Ref{Int64}         # Index of the best point so far in 'Y'.
    c     :: Ref{Float64}       # Constant of the model.
    g     :: AbstractVector     # Gradient of the model.
    gopt  :: AbstractVector     # Gradiente of the model in 'xopt'
    hq    :: AbstractVector     # Holds the explicit second derivatives of the quadratic model.
    pq    :: AbstractVector     # Holds parameters of the implicit second derivatives of the quadratic model.
    xbase :: AbstractVector     # Origin of the sample set.
    xopt  :: AbstractVector     # Best point so far (in terms of function values).
    fval  :: AbstractVector     # Set of the function values of the interpolation points.
    dst   :: AbstractVector     # Distances between 'xbase' and other interpolation points.
    Y     :: AbstractMatrix     # Set of interpolation points, shifted from the center of the sample set 'xbase'.
    BMAT  :: AbstractMatrix     # Holds the elements of 'Ξ', with the exception of its first column, 
                                # and the elements of 'Υ', with the exception of its first row and column.
    ZMAT  :: AbstractMatrix     # Holds the elements of 'Z', from the factorization 'Ω = ZZ^T'.

end