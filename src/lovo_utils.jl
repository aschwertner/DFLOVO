#-------------------------------------------------------------------------------
# Set of useful functions related to LOVO functions.
#-------------------------------------------------------------------------------

"""

    fmin_eval(func_list::Array{Function, 1}, r:: Int64, 
                y::Vector{Float64})

Computes the value of the objective function fmin(y) and an index belonging to
the set Imin(y). 

    - 'func_list': list containing the functions that determine the objective
    function fmin.

    - 'r': number of functions that make up the objective function fmin.

    - 'y': n-dimensional vector.
    
    Returns Returns the function value 'fmin_y' and the index 'imin_y'.

"""
function fmin_eval(
                    func_list::Array{Function, 1}, 
                    r:: Int64,
                    y::Vector{Float64}
                    )

    fmin_y = func_list[1](y)
    imin_y = 1
    for i = 2:r
        tmp = func_list[i](y)
        if tmp < fmin_y
            fmin_y = tmp
            imin_y = i
        end
    end

    return fmin_y, imin_y

end