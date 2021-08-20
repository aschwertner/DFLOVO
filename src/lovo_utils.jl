#-------------------------------------------------------------------------------
# Set of useful functions related to LOVO functions.
#-------------------------------------------------------------------------------

"""

    fmin_eval!(func_list::Array{Function, 1}, r:: Int64, 
                y::Vector{Float64}, imin_set::Vector{Bool})

Computes the value of the objective function fmin(y) and an index belonging to
the set I_{min}(y). 

    - 'func_list': list containing the functions that determine the objective
    function fmin.

    - 'r': number of functions that make up the objective function fmin.

    - 'y': n-dimensional vector.

The function modifies the argument:

    - 'imin_set': boolean vector with the indexes belonging to the I_{min}(y) set.
    
Returns the function value 'fmin_y' and the index 'imin_y'.

"""
function fmin_eval(
                    func_list::Array{Function, 1}, 
                    r::Int64,
                    y::Vector{Float64},
                    imin_set::Vector{Bool}
                    )

    fmin_y = func_list[1](y)
    imin_y = 1
    imin_set[1] = true 
    for i = 2:r
        tmp = func_list[i](y)
        if tmp < fmin_y
            fmin_y = tmp
            imin_y = i
            imin_set[1:(i - 1)] .= false
            imin_set[i] = true
        elseif tmp = fmin_y
            imin_set[i] = true
        else
            imin_set[i] = false
        end
    end

    return fmin_y, imin_y

end

"""

    fmin_partial_eval(func_list::Array{Function, 1}, r:: Int64, idx::Int64,
                        fi_y::Float64, y::Vector{Float64}, imin_set::Vector{Bool})

Computes the value of the objective function fmin(y) and an index belonging to
the set Imin(y), using the precalculated value of f_idx(y).

    - 'func_list': list containing the functions that determine the objective
    function fmin.

    - 'r': number of functions that make up the objective function fmin.

    - 'idx': index corresponding to the available function value.

    - 'fi_y': available function value.

    - 'y': n-dimensional vector.

The function modifies the argument:

    - 'imin_set': boolean vector with the indexes belonging to the I_{min}(y) set.
    
Returns the function value 'fmin_y' and the index 'imin_y'.

"""
function fmin_partial_eval(
                    func_list::Array{Function, 1}, 
                    r::Int64,
                    idx::Int64,
                    fi_y::Float64,
                    y::Vector{Float64},
                    imin_set::Vector{Bool}
                    )

    fmin_y = fi_y
    imin_y = idx
    imin_set[idx] = true
    for i = 1:r
        if i != idx
            tmp = func_list[i](y)
            if tmp < fmin_y
                fmin_y = tmp
                imin_y = i
                imin_set[1:(i - 1)] .= false
                imin_set[i] = true
            elseif tmp = fmin_y
                imin_set[i] = true
            else
                imin_set[i] = false
            end
        end
    end

    return fmin_y, imin_y

end

"""

    fi_eval(func_list::Array{Function, 1}, idx:: Int64, 
                y::Vector{Float64})

Computes the value of the function f_index(y).

    - 'func_list': list containing the functions that determine the objective
    function fmin.

    - 'idx': index of the funtion in 'func_list'.

    - 'y': n-dimensional vector.
    
Returns the function value.

"""
function fi_eval(
                    func_list::Array{Function, 1}, 
                    idx:: Int64,
                    y::Vector{Float64}
                    )

    return func_list[idx](y)

end