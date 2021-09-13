#-------------------------------------------------------------------------------
# Set of functions useful for the execution of the main algorithm. 
#-------------------------------------------------------------------------------

"""

    verify_initial_room(n::Int64, δ::Float64, a::Vector{Float64},
                        b::Vector{Float64})

Checks whether the bounds satisfy the conditions 'b' >= 'a' + 2 * 'δ'.

    - 'n': dimension of the search space.

    - 'δ': radius of the sample set.

    - 'a': n-dimensional vector with the lower bounds.

    - 'b': n-dimensional vector with the upper bounds.

Returns a Boolean value.

"""
function verify_initial_room(
                                n::Int64, 
                                δ::Float64, 
                                a::Vector{Float64}, 
                                b::Vector{Float64}
                                )

    for i = 1:n

        if ( b[i] - a[i] ) < ( 2.0 * δ )

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

The function modifies the argument:

    - 'x': n-dimensional vector (initial guess).

    - 'ao': n-dimensional vector with the shifted lower bounds.

    - 'bo': n-dimensional vector with the shifted upper bounds.

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

    - 'idx': index of the point in the sample set.

    - 'n': dimension of the search space.

    - 'a': n-dimensional vector with the lower bounds.

    - 'b': n-dimensional vector with the upper bounds.

    - 'ao': n-dimensional vector with the shifted lower bounds.

    - 'bo': n-dimensional vector with the shifted upper bounds.

    - 'xbase': n-dimensional vector (origin of the sample set).

    - 'Y': (n x m)-dimensional matrix (set of sample points).

The function modifies the argument:

    - 'x': n-dimensional vector.

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

function relative_reduction(
                            model::AbstractModel,
                            func_list::Array{Function, 1},
                            diff::Float64,
                            y::Vector{Float64}
                            )

    f_y = fi_eval(func_list, model.imin[], y)
    ρ = ( model.fval[model.kopt[] + 1] - f_y ) / ( diff )

    return ρ, f_y, model.imin[]

end

function relative_reduction!(
                            model::AbstractModel,
                            func_list::Array{Function, 1},
                            r::Int64,
                            diff::Float64,
                            y::Vector{Float64},
                            imin_set::Vector{Bool}
                            )

    fmin_y, idx_y = fmin_eval!(func_list, r, y, imin_set)
    ρ = ( model.fval[model.kopt[] + 1] - fmin_y ) / ( diff )

    return ρ, fmin_y, idx_y

end

"""

    print_warning(flag::Int64)

Prints information about the exit flag. 

    - 'flag': exit flag.

"""
function print_warning(
                        flag::Int64
                        )

    if flag == 1

        printstyled("Success: ", bold=true, color=:light_green)
        println("The algorithm succeeded.")

    elseif flag == -1

        printstyled("Warning: ", bold=true, color=:light_red)
        println("The maximum number of iterations has been reached.")

    elseif flag == -2

        printstyled("Warning: ", bold=true, color=:light_red)
        println("The maximum number of function evaluations has been reached.")

    elseif flag == -3

        printstyled("Warning: ", bold=true, color=:light_red)
        println("The calculated direction is not downhill for the model.")

    elseif flag == -11

        printstyled("Notice: ", bold=true, color=:light_yellow)
        println("A better value will be obtained when evaluating the function on 'xopt'.")

    else

        printstyled("Warning: ", bold=true, color=:light_red)
        println("Exit flag not specified.")

    end

end

"""

    print_iteration(countit::Int64, countf::Int64, flag::Int64, 
                    δ::Float64, Δ::Float64, fsave::Float64)
                            
Prints information about the iteration status.

    - 'countit': iteration count.

    - 'countf': function evaluation counter.

    - 'flag': iteration type flag.

    - 'δ': radius of the sample set.

    - 'Δ': radius of the trust-region.

    - 'fsave': least value of the objective function so far.
   
"""
function print_iteration(
                            full_calc::Bool,
                            flag::Int64,
                            countit::Int64,
                            countf::Int64,
                            imin_idx::Int64,
                            δ::Float64,
                            Δ::Float64,
                            fopt::Float64,
                            xopt::Vector{Float64}
                            )
        
    if flag == 1

        it_type = "criticality"

    elseif flag == 2

        it_type = "trust-region"

    elseif flag == 3

        it_type = "altmov"

    elseif flag == 4

        it_type = "altmov-cauchy"

    else

        it_type = "not especified"

    end

    if countit == 0

        println("--------------------------------------------------------------------------------")
        println("--------------------------------------------------------------------------------")

    end

    println("Iteration  : $(countit)")
    println("Func. eval.: $(countf)")
    println("δ          : $(δ)")
    println("Δ          : $(Δ)")
    println("It. type   : $(it_type)")
    println("I_min index: $(imin_idx)")
    println("Best point : $(xopt)")
    println("Func. val. : $(fopt)")
    println("Full ρ     : $(full_calc)")
    println("--------------------------------------------------------------------------------")

end

function print_info(
                    model::AbstractModel,
                    output::LOWDEROutput,
                    full_calc::Bool,
                    verbose::Int64,
                    exit_flag::Int64,
                    it_flag::Int64,
                    countit::Int64,
                    countf::Int64,
                    δ::Float64,
                    Δ::Float64
                    )
    
    if verbose != 0

        print_iteration( full_calc, it_flag, countit, countf, model.imin[], δ, Δ, model.fval[model.kopt[]], model.xopt)

        if verbose ≥ 2

            print_warning(exit_flag)

            if model.kopt[] != 0

                print_warning(-11)

            end

            if verbose == 3

                show_output(output)

            end

        end

    end

end

function projection_active_set!(
                                v::Vector{Float64}, 
                                active_set::Vector{Bool}, 
                                proj_v::Vector{Float64};
                                sym::Bool=false
                                )

    for i=1:length(v)

        if active_set[i]

            proj_v[i] = 0.0

        else

            if sym

                proj_v[i] = - v[i]

            else

                proj_v[i] = v[i]

            end

        end
        
    end

end