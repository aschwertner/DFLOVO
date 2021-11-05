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

Returns a boolean value.

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

    relative_reduction(model::AbstractModel, func_list::Array{Function, 1},
                        pred_red::Float64, y::Vector{Float64})

Computes the simplified relative reduction (using f_{i} instead of fmin) in the point 'y'.

    - 'model': model of LinearModel or QuadraticModel type.

    - 'func_list': list containing the functions that determine the objective
    function fmin.

    - 'pred_red': precalculated predicted reduction.

    - 'y': n-dimensional vector (point of interest).

Returns the value of the simplified relative reduction, the real reduction, the value of f_i(y) and the index 'model.imin'.

"""
function relative_reduction(
                            model::AbstractModel,
                            func_list::Array{Function, 1},
                            pred_red::Float64,
                            y::Vector{Float64}
                            )

    f_y = fi_eval(func_list, model.imin[], y)
    real_red = model.fval[model.kopt[] + 1] - f_y
    ρ = real_red / pred_red

    return ρ, real_red, f_y, model.imin[]

end

"""

    relative_reduction(model::AbstractModel, func_list::Array{Function, 1}, r::Int64,
                        pred_red::Float64, y::Vector{Float64}, imin_set::Vector{Bool})

Computes the relative reduction obtained in the point 'y'.

    - 'model': model of LinearModel or QuadraticModel type.

    - 'func_list': list containing the functions that determine the objective
    function fmin.

    - 'r': number of functions that make up the objective function fmin.

    - 'pred_red': precalculated predicted reduction.

    - 'y': n-dimensional vector (point of interest).

The function modifies the argument:

    - 'imin_set': boolean vector with the indexes belonging to the I_{min}(y) set.

Returns the value of the relative reduction, the real reduction, the value of fmin(y) and an index 'idx_y' belonging to Imin(y).

"""
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
    ρ = real_red / pred_red

    return ρ, real_red, fmin_y, idx_y

end

"""

    reconstruct_original_point!(model::AbstractModel, idx::Int64, a::Vector{Float64}, b::Vector{Float64},
                                ao::Vector{Float64}, bo::Vector{Float64}, x::Vector{Float64})
    
Constructs the model based in the given information.

    - 'model': model of LinearModel or QuadraticModel type.

    - 'idx': position of the point in 'Y' set.

    - 'a': n-dimensional vector with the lower bounds.

    - 'b': n-dimensional vector with the upper bounds.

    - 'ao': n-dimensional vector with the shifted lower bounds.

    - 'bo': n-dimensional vector with the shifted upper bounds.

The function modifies the argument:

    - 'x': n-dimensional vector (point of interest).    

"""
function reconstruct_original_point!(
                                    model::AbstractModel,
                                    idx::Int64,
                                    a::Vector{Float64},
                                    b::Vector{Float64},
                                    ao::Vector{Float64},
                                    bo::Vector{Float64},
                                    x::Vector{Float64}
                                    )

    for i = 1:model.n
     
        if model.Y[idx, i] == ao[i]

            x[i] = a[i]

        elseif model.Y[idx, i] == bo[i]

            x[i] = b[i]

        else

            x[i] = min( max( a[i], model.xbase[i] + model.Y[idx, i] ), b[i] )

        end

    end

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

    print_iteration(it_flag::Symbol, status_flag::Symbol, full_calc::Bool, nit::Int64,
                        nf::Int64, imin_idx::Int64, δ::Float64, Δ::Float64, π::Float64,
                        ρ::Float64, pred_red::Float64, real_red::Float64, fopt::Float64,
                        xopt::Vector{Float64}, d::Vector{Float64})

Prints information about the iterations.

    - 'it_flag': iteration flag.

    - 'status_flag': status flag.

    - 'full_calc': indicates if ρ was full calculated.

    - 'nit': number of iterations.

    - 'nf': number of function evaluations.

    - 'imin_idx': index corresponding to the available function value.

    - 'δ': sample set radius.

    - 'Δ': trust-region radius.

    - 'π': stationarity measure.
 
    - 'ρ': relative reduction.

    - 'pred_red': predicted reduction.

    - 'real_red': real reduction.

    - 'fopt': function value in 'xopt'

    - 'xopt': best point so far.

    - 'd': last computed direction.

"""
function print_iteration(
                            it_flag::Symbol,
                            status_flag::Symbol,
                            full_calc::Bool,
                            nit::Int64,
                            nf::Int64,
                            imin_idx::Int64,
                            δ::Float64,
                            Δ::Float64,
                            π::Float64,
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

    if ( nit == 0 && it_flag == :max_evaluations )

        println("Iteration     : $( nit )")
        println("Func. eval.   : $( nf )")
        println("Radius δ      : $( δ )")
        println("Radius Δ      : $( Δ )")
        println("I_min index   : $( imin_idx )")
        println("Best point    : $( xopt )")
        println("Func. val.    : $( fopt )")
        println("Iter. type    : $( it_flag )")
        println("--------------------------------------------------------------------------------")

    elseif it_flag == :nonspecified

        println("Iteration     : $( nit )")
        println("Func. eval.   : $( nf )")
        println("Radius δ      : $( δ )")
        println("Radius Δ      : $( Δ )")
        println("Stationarity π: $( π )")
        println("I_min index   : $( imin_idx )")
        println("Best point    : $( xopt )")
        println("Func. val.    : $( fopt )")
        println("Iter. type    : $( it_flag )")
        println("--------------------------------------------------------------------------------")

    elseif it_flag == :criticality

        println("Iteration     : $( nit )")
        println("Func. eval.   : $( nf )")
        println("Radius δ      : $( δ )")
        println("Radius Δ      : $( Δ )")
        println("Stationarity π: $( π )")
        println("I_min index   : $( imin_idx )")
        println("Best point    : $( xopt )")
        println("Func. val.    : $( fopt )")
        println("Iter. type    : $( it_flag )")
        println("Step status   : $( status_flag )")
        println("Direction d   : $( d )")
        println("--------------------------------------------------------------------------------")

    else

        println("Iteration     : $( nit )")
        println("Func. eval.   : $( nf )")
        println("Radius δ      : $( δ )")
        println("Radius Δ      : $( Δ )")
        println("Stationarity π: $( π )")
        println("I_min index   : $( imin_idx )")
        println("Best point    : $( xopt )")
        println("Func. val.    : $( fopt )")
        println("Iter. type    : $( it_flag )")
        println("Step status   : $( status_flag )")
        println("Direction d   : $( d )")
        println("Full ρ        : $( full_calc )")
        println("Rel. reduc. ρ : $( ρ )")
        println("Pred. reduc.  : $( pred_red )")
        println("Real reduc.   : $( real_red )")
        println("--------------------------------------------------------------------------------")

    end

end

"""

    print_info(model::AbstractModel, output::LOWDEROutput, exit_flag::Symbol, it_flag::Symbol,
                status_flag::Symbol, verbose::Int64, nit::Int64, nf::Int64, δ::Float64, Δ::Float64,
                π::Float64, full_calc::Bool, pred_red::Float64, real_red::Float64, ρ::Float64,
                d::Vector{Float64})

Prints information about the last iteration of LOWDER.

    - 'model': model of LinearModel or QuadraticModel type.

    - 'output': object of LOWDEROutput type.

    - 'it_flag': iteration flag.

    - 'status_flag': status flag.

    - 'full_calc': indicates if fmin was calculated at 'xopt'.

    - 'verbose': level of vebosity.

    - 'δ': sample set radius.

    - 'Δ': trust-region radius.

    - 'π': stationarity measure.

    - 'ρ': relative reduction.
 
    - 'pred_red': predicted reduction.

    - 'real_red': real reduction.

    - 'd': last computed direction.

"""
function print_info(
                    model::AbstractModel,
                    output::LOWDEROutput,
                    it_flag::Symbol,
                    status_flag::Symbol,
                    full_calc::Bool,
                    verbose::Int64,
                    δ::Float64,
                    Δ::Float64,
                    π::Float64,
                    ρ::Float64,
                    pred_red::Float64,
                    real_red::Float64,
                    d::Vector{Float64}
                    )
    
    if verbose != 0

        print_iteration(it_flag, status_flag, full_calc, output.iter, output.nf, output.index, δ, Δ, π, ρ, pred_red, real_red, output.f, output.solution, d)

        if verbose ≥ 2

            print_warning(output.status)

            if !( model.kopt[] == 0 || full_calc ) 

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

function step_projection!(
                            model::AbstractModel,
                            a::Vector{Float64},
                            b::Vector{Float64},
                            x::Vector{Float64},
                            d::Vector{Float64}
                            )

    for i = 1:model.n

        x[i] = max( a[i], min( model.xopt[i] + d[i], b[i] ) )
        d[i] = x[i] - model.xopt[i]

    end

end

function compute_alpha_linear(
                                Δ::Float64,
                                a::Vector{Float64},
                                b::Vector{Float64},
                                x::Vector{Float64},
                                d::Vector{Float64},
                                s::Vector{Float64}
                                )

    roots = []
    α_j = Inf
    α_B = Inf
    idx = 0

    # Computes α_B, the largest number such that 'a' ≤ 'x' + 'd' + 'α'*'s' ≤ 'b'.
    for j=1:length(x)

        if s[j] < 0.0

            α_j = ( a[j] - x[j] - d[j] ) / s[j]

        elseif s[j] > 0.0

            α_j = ( b[j] - x[j] - d[j] ) / s[j]

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

function new_search_direction!(
                                pdTpd::Float64,
                                pdTpg::Float64,
                                pgTpg::Float64,
                                proj_d::Vector{Float64},
                                proj_grad::Vector{Float64}
                                )

    roots = []
    α = 0.0
    β = 0.0

    if pdTpg ≈ 0.0

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

        if ( α == NaN ) || ( α * pdTpd ≈ 0.0 )

            @. proj_grad = 0.0

        else

            β = - pgTpg / ( α * pdTpd )

            if ( α * pdTpg + β * pgTpg ) < 0.0

                @. proj_grad = α * proj_d + β * proj_grad

            else

                α = minimum(roots)

                if ( α == NaN ) || ( α * pdTpd ≈ 0.0 )

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

function solve_quadratic!(
                            a::Float64,
                            b::Float64,
                            c::Float64,
                            roots::Vector{Any}
                            )

    if a ≈ 0.0

        root = - c / b

        if !( isinf(root) )

            append!(roots, root)

        elseif ( b ≈ 0.0 ) && ( c ≈ 0.0 )

            append!(roots, 0.0)

        else

            append!(roots, NaN)

        end

    else

        b = b / a
        c = c / a

        if c ≈ 0.0

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

function binary_search(
                        lower_value::Float64,
                        upper_value::Float64,
                        stop_condition::Function,
                        ε::Float64
                        )

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

function cond_θB(
                    θ::Float64,
                    a::Vector{Float64}, 
                    b::Vector{Float64},
                    xopt::Vector{Float64},
                    d::Vector{Float64},
                    proj_d::Vector{Float64},
                    s::Vector{Float64}
                    )

    for i=1:length(a)

        xnewi = xopt[i] + d[i] - proj_d[i] + cos(θ) * proj_d[i] + sin(θ) * s[i]

        if ( xnewi < a[i] ) || ( xnewi > b[i] )

            return false

        end

    end

    return true

end

function cond_θQ_linear(
                        θ::Float64, 
                        gTd::Float64, 
                        gTpd::Float64, 
                        gTs::Float64
                        )

    if ( gTd - gTpd + cos(θ) * gTpd + sin(θ) * gTs ) < 0.0

        return true

    else

        return false

    end

end

function update_active_set!(
                            a::Vector{Float64},
                            b::Vector{Float64},
                            xopt::Vector{Float64},
                            d::Vector{Float64},
                            active_set::Vector{Bool}
                            )

    for i=1:length(a)

        if !( active_set[i] )
            
            xnewi = xopt[i] + d[i]

            if ( a[i] == xnewi ) || ( xnewi == b[i] )

                active_set[i] == true

            end

        end

    end

end

function save_info!(
                    model::AbstractModel,
                    it_flag::Symbol,
                    status_flag::Symbol,
                    full_calc::Bool,
                    nit::Int64,
                    nf::Int64,
                    δ::Float64,
                    Δ::Float64,
                    π::Float64,
                    ρ::Float64,
                    pred_red::Float64,
                    real_red::Float64,
                    d::Vector{Float64},
                    file::IO
                    )
    
    text = "$( nit ) $( nf ) $( model.imin[] ) $( model.kopt[] ) $( δ ) $( Δ ) $( π ) $( it_flag ) $( status_flag ) $( full_calc ) $( ρ ) $( pred_red ) $( real_red ) $( model.c[] )"

    print( file, text )

    for i = 1:model.n

        print( file, " $( model.g[i] )" )

    end

    for i = 1:model.n

        print( file, " $( model.xbase[i] )" )

    end

    for i = 1:model.n

        print( file, " $( model.xopt[i] )" )

    end

    for i = 1:model.m

        print( file, " $( model.fval[i] )" )

    end

    for i = 1:(model.m - 1)

        for j = 1:model.n

            print( file, " $( model.Y[i, j] )" )

        end

    end

    for i = 1:(model.n - 1)

        print( file, " $( d[i] )" )

    end

    println( file, " $( d[model.n] )" )
    
end

function max_distance(
                        model::AbstractModel
                        )

    max_dst = - Inf

    for i = 1:( model.m - 1 )

        aux = norm( model.Y[i, :] )

        if max_dst ≤ aux

            max_dst = aux

        end

    end

    return max_dst

end

function choose_index_altmov(
                                model::AbstractModel
                            )
    
    ( val, idx_t ) = findmax( model.dst )

    if idx_t == model.kopt[]

        idx_t = 0

    end

    return idx_t

end