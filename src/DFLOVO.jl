# Turn off precompilation during development
__precompile__(false)

# Main file implementing DFLOVO in Julia

module DFLOVO

    import Base: (*)
    using LinearAlgebra

    #Adds files
    include("lovo_utils.jl")
    include("utils.jl")

end
