#-------------------------------------------------------------------------------
# Set of functions useful for the execution of the main algorithm. 
#-------------------------------------------------------------------------------

"""

    verify_initial_room(n::Int64, Δ::Float64, a::Vector{Float64},
                        b::Vector{Float64})

Checks whether the bounds satisfy the conditions 'b' >= 'a' + 2 * 'Δ'.

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

"""

    correct_guess_bounds!(n::Int64, δ::Float64, a::Vector{Float64},
                            b::Vector{Float64}, x::Vector{Float64},
                            ao::Vector{Float64}, bo::Vector{Float64})

Fixes initial guess 'x' and sets the lower and upper bound vectors 'ao' and 'bo'
shifted to the origin 'x'.

    - 'n': dimension of the search space.

    - 'δ': radius of the sample set.

    - 'a': n-dimensional vector with the lower bounds.

    - 'b': n-dimensional vector with the upper bounds.

    - 'x': n-dimensional vector (initial guess).

    - 'ao': n-dimensional vector with the shifted lower bounds.

    - 'bo': n-dimensional vector with the shifted upper bounds.

Modifies the vectors 'x', 'ao', 'bo'.

"""
function correct_guess_bounds!( 
                                n::Int64, 
                                δ::Float64, 
                                a::Vector{Float64}, 
                                b::Vector{Float64},
                                x::Vector{Float64},
                                ao::Vector{Float64}, 
                                bo::Vector{Float64}
                                )

    for i=1:n
        ao[i] = a[i] - x[i]
        bo[i] = b[i] - x[i]
        if ao[i] >= - δ
            if ao[i] >= 0.0
                x[i] = a[i]
                ao[i] = 0.0
                bo[i] = b[i] - a[i]
            else
                x[i] = a[i] + δ
                ao[i] = - δ
                bo[i] = max( b[i] - x[i], δ )
            end
        elseif bo[i] <= δ
            if bo[i] <= 0.0
                x[i] = b[i]
                ao[i] = a[i] - b[i]
                bo[i] = 0.0
            else
                x[i] = b[i] - δ
                ao[i] = min( a[i] - x[i], - δ )
                bo[i] = δ
            end
        end
    end

end

"""

    reconstruct_original_point!(idx::Int64, n::Int64, a::Vector{Float64}, 
                                b::Vector{Float64}, ao::Vector{Float64}, 
                                bo::Vector{Float64}, xbase::Vector{Float64},
                                Y::Matrix{Float64}, x::Vector{Float64})

Reconstructs the original point given its position 'idx' in the sample set  'Y'.

    - 'idx': index of the funtion in 'func_list'.

    - 'n': dimension of the search space.

    - 'a': n-dimensional vector with the lower bounds.

    - 'b': n-dimensional vector with the upper bounds.

    - 'ao': n-dimensional vector with the shifted lower bounds.

    - 'bo': n-dimensional vector with the shifted upper bounds.

    - 'xbase': n-dimensional vector (origin of the sample set).

    - 'Y': (n x m)-dimensional matrix (set of sample points).

    - 'x': n-dimensional vector.

Modifies the vector 'x' to become the desired point.

"""
function reconstruct_original_point!(
                                        idx::Int64,
                                        n::Int64,
                                        a::Vector{Float64}, 
                                        b::Vector{Float64},
                                        ao::Vector{Float64}, 
                                        bo::Vector{Float64},
                                        xbase::Vector{Float64},
                                        Y::Matrix{Float64},
                                        x::Vector{Float64}
                                        )
    
    for i=1:n
        x[i] = min( max( a[i], xbase[i] + Y[i, idx] ), b[i] )
        if Y[i, idx] == ao[i]
            x[i] = a[i]
        elseif Y[i, idx] == bo[i]
            x[i] = b[i]
        end
    end

end