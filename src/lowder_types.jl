#-------------------------------------------------------------------------------
# LOWDER types
#-------------------------------------------------------------------------------

export AbstractModel, LOWDEROutput

abstract type AbstractModel end

struct LOWDEROutput
   
    status     :: Symbol
    true_val   :: Bool
    iter       :: Int64
    nf         :: Int64
    index      :: Int64
    f          :: Float64
    solution   :: Vector{Float64}

end
