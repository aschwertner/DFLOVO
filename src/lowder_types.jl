#-------------------------------------------------------------------------------
# LOWDER types
#-------------------------------------------------------------------------------

export AbstractModel, LOWDEROutput

abstract type AbstractModel end

"""

Define the type for LOWDER output.

    - 'status': solution status.
    - 'true_val': indicates whether the value 'f' is the true value of the function 'fmin' for the point 'solution'.
    - 'iter': number of iterations.
    - 'nf': number of function evaluations.
    - 'index': index of the function in 'func_list'.
    - 'f': objective function value.
    - 'solution': solution found by LOWDER.

"""
struct LOWDEROutput
   
    status     :: Symbol
    true_val   :: Bool
    iter       :: Int64
    nf         :: Int64
    index      :: Int64
    f          :: Float64
    solution   :: Vector{Float64}

end
