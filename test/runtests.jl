#-------------------------------------------------------------------------------
# Run test files
#-------------------------------------------------------------------------------

# Adds modules
using LOWDER
using Test
using LinearAlgebra

# Adds files
include("test_lovo_utils.jl")
include("test_utils.jl")

include("test_linear_models.jl")
