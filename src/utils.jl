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
                            pred_red::Float64,
                            y::Vector{Float64}
                            )

    f_y = fi_eval(func_list, model.imin[], y)
    real_red = model.fval[model.kopt[] + 1] - f_y
    ρ = ( real_red ) / ( pred_red )

    return ρ, real_red, f_y, model.imin[]

end

function relative_reduction!(
                            model::AbstractModel,
                            func_list::Array{Function, 1},
                            r::Int64,
                            pred_red::Float64,
                            y::Vector{Float64},
                            imin_set::Vector{Bool}
                            )

    fmin_y, idx_y = fmin_eval!(func_list, r, y, imin_set)
    real_red = model.fval[model.kopt[] + 1] - fmin_y
    ρ = ( real_red ) / ( pred_red )

    return ρ, real_red, fmin_y, idx_y

end

"""

    print_warning(flag::Int64)

Prints information about the exit flag. 

    - 'flag': exit flag.

"""
function print_warning(
                        exit_flag::Symbol
                        )

    if exit_flag == :success

        printstyled("Success: ", bold=true, color=:light_green)
        println("The algorithm succeeded.")

    elseif exit_flag == :max_iterations

        printstyled("Warning: ", bold=true, color=:light_red)
        println("The maximum number of iterations has been reached.")

    elseif exit_flag == :max_evaluations

        printstyled("Warning: ", bold=true, color=:light_red)
        println("The maximum number of function evaluations has been reached.")

    elseif exit_flag == :nondescent

        printstyled("Warning: ", bold=true, color=:light_red)
        println("The calculated direction is not descent for the model.")

    elseif exit_flag == :better_point

        printstyled("Notice: ", bold=true, color=:light_yellow)
        println("A better value will be obtained when evaluating the function on 'xopt'.")

    else

        printstyled("Warning: ", bold=true, color=:light_red)
        println("Exit flag not specified.")

    end
    println("--------------------------------------------------------------------------------")

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
                            it_flag::Symbol,
                            nit::Int64,
                            nf::Int64,
                            imin_idx::Int64,
                            δ::Float64,
                            Δ::Float64,
                            fopt::Float64,
                            xopt::Vector{Float64}
                            )
        
    if nit == 0

        println("--------------------------------------------------------------------------------")
        println("--------------------------------------------------------------------------------")

    end

    println("Iteration    : $(nit)")
    println("Func. eval.  : $(nf)")
    println("δ            : $(δ)")
    println("Δ            : $(Δ)")
    println("I_min index  : $(imin_idx)")
    println("Best point   : $(xopt)")
    println("Func. val.   : $(fopt)")
    println("Iter. type   : $(it_flag)")
    println("--------------------------------------------------------------------------------")

end

function print_iteration(
                            it_flag::Symbol,
                            full_calc::Bool,
                            nit::Int64,
                            nf::Int64,
                            imin_idx::Int64,
                            δ::Float64,
                            Δ::Float64,
                            ρ::Float64,
                            pred_red::Float64,
                            real_red::Float64,
                            fopt::Float64,
                            xopt::Vector{Float64},
                            d::Vector{Float64}
                            )
        
    if nit == 0

        println("--------------------------------------------------------------------------------")
        println("--------------------------------------------------------------------------------")

    end

    println("Iteration    : $(nit)")
    println("Func. eval.  : $(nf)")
    println("δ            : $(δ)")
    println("Δ            : $(Δ)")
    println("I_min index  : $(imin_idx)")
    println("Best point   : $(xopt)")
    println("Func. val.   : $(fopt)")
    println("Iter. type   : $(it_flag)")
    println("Direction d  : $(d)")
    println("Full ρ       : $(full_calc)")
    println("Rel. reduc. ρ: $(ρ)")
    println("Pred. reduc. : $(pred_red)")
    println("Real reduc.  : $(real_red)")
    println("--------------------------------------------------------------------------------")

end

function print_info(
                    model::AbstractModel,
                    output::LOWDEROutput,
                    exit_flag::Symbol,
                    it_flag::Symbol,
                    verbose::Int64,
                    nit::Int64,
                    nf::Int64,
                    δ::Float64,
                    Δ::Float64,
                    full_calc::Bool,
                    pred_red::Float64,
                    real_red::Float64,
                    ρ::Float64,
                    d::Vector{Float64}
                    )
    
    if verbose != 0

        if it_flag == :nonspecified

            print_iteration(it_flag, nit, nf, model.imin[], δ, Δ, model.fval[model.kopt[] + 1], model.xopt)
        
        else

            print_iteration(it_flag, full_calc, nit, nf, model.imin[], δ, Δ, ρ, pred_red, real_red, model.fval[model.kopt[] + 1], model.xopt, d)

        end

        if verbose ≥ 2

            print_warning(exit_flag)

            if model.kopt[] != 0

                print_warning(:better_point)

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

function compute_alpha_linear(Δ, a, b, x, d, s)

    roots = []
    α_j = Inf
    α_B = Inf
    idx = 0

    # Computes α_B, the largest number such that 'a' ≤ 'x' + 'd' + 'α'*'s' ≤ 'b'.
    for j=1:length(x)

        if s[i] < 0.0

            α_j = ( a[i] - x[i] - d[i] ) / s[i]

        elseif s[i] > 0.0

            α_j = ( b[i] - x[i] - d[i] ) / s[i]

        end

        if α_j < α_B

            α_B = α_j
            idx = j

        end

    end

    # Computes α_Δ, the largest number such that ||'d' + 'α'*'s'|| ≤ 'Δ'.

    a_2 = dot(s, s)
    a_1 = 2.0 * dot(d, s)
    a_0 = dot(d, d) - Δ ^ 2.0

    solve_quadratic!(a_2, a_1, a_0, roots)

    α_Δ = maximum(roots)

    if ( α_Δ == NaN ) || ( α_Δ < 0.0 )

        return α_B, 2, idx

    else

        α = min(α_B, α_Δ)

        if α == α_Δ 

            return α, true, 0

        else

            return α, false, idx

        end

    end

end

function new_search_direction!(pdTpd, pdTpg, pgTpg, proj_d, proj_grad)

    roots = []
    α = 0.0
    β = 0.0

    if iszero( pdTpg )

        # If P_{I}(d)'P_{I}(∇model(xk+d)) = 0, then the new direction is a multiple of P_{I}(∇Q(xk+d)).

        solve_quadratic!( pgTpg, 0.0, pdTpd, roots)

        β = minimum(roots)

        if ( β == NaN ) || ( β > 0.0 )

            @. proj_grad = 0.0
        
        else

            @. proj_grad *= β

        end

    else

        # Otherwise, the new search direction is a linear combination of P_{I}(d) and P_{I}(∇model(xk+d)).

        aux = ( pdTpd ^ 2.0 * pgTpg ) / pdTpg ^ 2.0 - pdTpd

        solve_quadratic!( aux, 0.0, - pdTpd, roots )

        α = maximum(roots)

        if ( α == NaN ) || ( iszero( α * pdTpd ) )

            @. proj_grad = 0.0

        else

            β = - pgTpg / ( α * pdTpd )

            if ( α * pdTpg + β * pgTpg ) < 0.0

                @. proj_grad = α * proj_d + β * proj_grad

            else

                α = minimum(roots)

                if ( α == NaN ) || ( iszero( α * pdTpd ) )

                    @. proj_grad = 0.0
        
                else

                    β = - pgTpg / ( α * pdTpd )

                    if ( α * pdTpg + β * pgTpg ) < 0.0

                        @. proj_grad = α * proj_d + β * proj_grad
        
                    else

                        @. proj_grad = 0.0

                    end

                end

            end

        end

    end

end

function solve_quadratic!(a, b, c, roots)

    if iszero(a)

        root = - c / b

        if !( isinf(root) )

            append!(roots, root)

        elseif ( iszero(b) ) && ( iszero(c) )

            append!(roots, 0.0)

        else

            append!(roots, NaN)

        end

    else

        b = b / a
        c = c / a

        if iszero(c)

            append!(roots, 0.0)
            append!(roots, -b)

        else

            Δ = b ^ 2.0 - 4.0 * c

            if isinf(Δ)

                root_1 = -b
                root_2 = c / root_1

                if isinf(root_2)

                    append!(roots, root_1)

                else

                    append!(roots, root_1)
                    append!(roots, root_2)
                
                end

            else

                if Δ < 0.0

                    append!(roots, NaN)

                elseif Δ == 0.0

                    append!(roots, - 0.5 * b)

                else
                    
                    if b < 0.0
                    
                        root_1 = 0.5 * ( -b + sqrt(Δ) )

                    else

                        root_1 = 0.5 * ( - b - sqrt(Δ) )

                    end

                    root_2 = c / root_1

                    if isinf(root_2)

                        append!(roots, root_1)

                    else

                        append!(roots, root_1)
                        append!(roots, root_2)
                
                    end

                end

            end

        end

    end

end

function binary_search(lower_value, upper_value, stop_condition, ε)

    if !( stop_condition(lower_value) )

        return NaN

    elseif stop_condition(upper_value)

        return upper_value

    else

        while true

            mean_value = (lower_value + upper_value) / 2.0

            if stop_condition(mean_value)

                if (upper_value - mean_value) ≤ ε

                    return mean_value

                else

                    lower_value = mean_value

                end

            else

                if (mean_value - lower_value) ≤ ε

                    return lower_value

                else

                    upper_value = mean_value

                end
                
            end

        end

    end

end

function cond_θB(θ, a, b, xopt, d, proj_d, s)

    for i=1:length(a)

        xnewi = xopt[i] + d[i] - proj_d[i] + cos(θ) * proj_d[i] + sin(θ) * s[i]

        if ( xnewi < a[i] ) || ( xnewi > b[i] )

            return false

        end

    end

    return true

end