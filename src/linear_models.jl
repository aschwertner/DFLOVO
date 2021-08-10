#-------------------------------------------------------------------------------
# Set of useful functions related to LinearModel struct.
#-------------------------------------------------------------------------------

struct LinearModel
   
    imin  :: Int64
    kopt  :: Int64
    fbase :: Float64
    c     :: Float64
    g     :: AbstractVector
    xbase :: AbstractVector
    xopt  :: AbstractVector
    fval  :: AbstractVector    
    dst   :: AbstractVector

end

