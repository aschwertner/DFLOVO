#-------------------------------------------------------------------------------
# LOWDER types
#-------------------------------------------------------------------------------

export AbstractModel, LOWDEROutput

abstract type AbstractModel end

struct LOWDEROutput
   
    iter       :: Int64
    index      :: Int64
    nf         :: Int64
    status     :: Int64
    true_val   :: Bool
    f          :: Float64
    solution   :: Vector{Float64}

end
