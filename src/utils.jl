#-------------------------------------------------------------------------------
# Set of functions useful for the execution of the main algorithm. 
#-------------------------------------------------------------------------------

"""

    verify_initial_room(n::Int64, Δ::Float64, a::Vector{Float64},
                        b::Vector{Float64})

Checks whether the limits satisfy the conditions 'b' >= 'a' + 2 * 'Δ'.

    - 'n': dimension of the search space.

    - 'Δ': radius of the trust-region.

    - 'a': n-dimensional vector with the lower bounds.

    - 'b': n-dimensional vector with the upper bounds.

Returns a Boolean value.

"""
function verify_initial_room(
                                n::Int64, 
                                Δ::Float64, 
                                a::Vector{Float64}, 
                                b::Vector{Float64}
                                )

    for i = 1:n
        if ( b[i] - a[i] ) < ( 2.0 * Δ )
            return false 
        end
    end

    return true

end