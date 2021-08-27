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

"""

    construct_initial_set!(func_list::Array{Function,1}, n::Int64, m::Int64,
                            imin::Int64, maxfun::Int64, fbase::Float64,
                            δ::Float64, a::Vector{Float64}, b::Vector{Float64},
                            ao::Vector{Float64}, bo::Vector{Float64},
                            xbase::Vector{Float64}, countf::Int64,
                            xopt::Vector{Float64}, fval::Vector{Float64},
                            gopt::Vector{Float64}, hq::Vector{Float64},
                            BMAT::Matrix{Float64}, ZMAT::Matrix{Float64},
                            Y::Matrix{Float64}, x::Vector{Float64})
                            
Builds sample set 'Y', vector of function values 'fval', optimal point 'xopt'
and the first interpolation model.

    - 'func_list': list containing the functions that determine the objective
    function fmin.

    - 'n': dimension of the search space.

    - 'm': number of interpolation conditions (points).

    - 'imin': index belonging to the set Imin(xbase).

    - 'maxfun': maximum number of function evaluations allowed.

    - 'fbase': objective function value in 'xbase'.

    - 'δ': radius of the sample set.

    - 'a': n-dimensional vector with the lower bounds.

    - 'b': n-dimensional vector with the upper bounds.

    - 'ao': n-dimensional vector with the shifted lower bounds.

    - 'bo': n-dimensional vector with the shifted upper bounds.

    - 'xbase': n-dimensional vector (origin of the sample set).

The function modifies the arguments:

    - 'countf': function evaluation counter.

    - 'xopt': n-dimensional vector (optimal point so far).

    - 'fval': m-dimensional vector (function values of the interpolation
    points).

    - 'gopt': n-dimensional vector (gradient of the quadratic model at 
    'xbase + xopt').

    - 'hq': (n(n+1)/2)-dimensional vector (explicit second derivatives of the
    quadratic model).

    - 'BMAT': ((m + n) x n)-dimensional matrix (elements of Ξ and Υ).

    - 'ZMAT': (m x (m - n - 1))-dimensional matrix (elements of Z).

    - 'Y': (n x m)-dimensional matrix (set of sample points).

    - 'x': n-dimensional vector.

Returns the index of 'xopt' in the set 'Y' and 'countf'.

"""
function construct_initial_set!(
                                func_list::Array{Function,1},
                                n::Int64,
                                m::Int64,
                                imin::Int64,
                                maxfun::Int64,
                                fbase::Float64,
                                δ::Float64,
                                a::Vector{Float64}, 
                                b::Vector{Float64},
                                ao::Vector{Float64}, 
                                bo::Vector{Float64},
                                xbase::Vector{Float64},
                                countf::Int64,
                                xopt::Vector{Float64},
                                fval::Vector{Float64},
                                gopt::Vector{Float64},
                                hq::Vector{Float64},
                                BMAT::Matrix{Float64},
                                ZMAT::Matrix{Float64},
                                Y::Matrix{Float64},
                                x::Vector{Float64}
                                )
    
    fval[1] = fbase
    kopt = 1
    aux = 1.0 / δ ^ 2.0
    new_countf = countf

    for i = 0:min( m - 1, maxfun - countf )
        j = i - n
        k = i + 1
        if i <= (2 * n)
            if (i >= 1) && ( i <= n )
                α = δ
                if bo[i] == 0.0
                    α *= - 1.0
                end
                Y[i, k] = α
            elseif i > n
                α = Y[j, k - n]
                β = - δ
                if ao[j] == 0.0
                    β = min(2.0 * δ, bo[j])
                elseif bo[j] == 0.0
                    β = max(- 2.0 * δ, ao[j])
                end
                Y[j, k] = β
            end
        else
            tmp = div(i - n - 1, n)
            p_j = i - tmp * n - n
            q_j = p_j + tmp
            if q_j > n
                tmp = p_j
                p_j = q_j - n
                q_j = tmp
            end
            Y[p_j, k] = Y[p_j, p_j + 1]
            Y[q_j, k] = Y[q_j, q_j + 1]
        end

        if i > 0
            reconstruct_original_point!(k, n, a, b, ao, bo, xbase, Y, x)
            fk = fi_eval(func_list, imin, x)
            new_countf += 1
        end

        if fk < fval[kopt]
            kopt = k
        end

        if k <= ( 2 * n + 1 )
            if ( k >= 2 ) && ( k <= n + 1 )
                gopt[i] = (fk - fbase) / α
                if m < ( k + n )
                    BMAT[1, i] = - 1.0 / α
                    BMAT[k, i] = 1.0 / α
                    BMAT[m + i, i] = - 0.5 * δ ^ 2.0
                end
            elseif k >= n + 2
                idx_h = convert(Int64, j * ( j + 1 ) / 2)
                tmp = (fk - fbase) / β
                diff = β - α
                hq[idx_h] = 2.0 * ( tmp - gopt[j] ) / diff
                gopt[j] = ( gopt[j] * β - tmp * α ) / diff
                if α * β < 0.0
                    if fk < fval[k - n]
                        fval[k] = fval[k - n]
                        fval[k - n] = fk
                        if kopt == k
                            kopt = k - n
                        end
                        Y[j, k - n] = β
                        Y[j, k] = α
                    end
                end
                BMAT[1, j] = - ( α + β ) / ( α * β )
                BMAT[k, j] = - 0.5 / Y[j, k - n]
                BMAT[k - n, j] = - BMAT[1, j] - BMAT[k, j]
                ZMAT[1, j] = sqrt( 2.0 ) / ( α * β )
                ZMAT[k, j] = sqrt( 0.5 ) / δ ^ 2.0
                ZMAT[k - n, j] = - ZMAT[1, j] - ZMAT[k, j]
            end
        else
            idx_h = convert(Int64, q_j * ( q_j - 1 ) / 2) + p_j
            ZMAT[1, j] = aux
            ZMAT[k, j] = aux
            ZMAT[q_j + 1, j] = - aux
            ZMAT[p_j + 1, j] = - aux
            tmp = Y[q_j, k] * Y[p_j, k]
            hq[idx_h] = ( fbase - fval[q_j + 1] - fval[p_j + 1] + fk ) / tmp
        end
    end

    reconstruct_original_point!(kopt, n, a, b, ao, bo, xbase, Y, xopt)

    return kopt, new_countf

end

"""

    update_gopt!(n::Int64, m::Int64, r::Int64, countf::Int64,
                    xopt::Vector{Float64}, hq::Vector{Float64},
                    pq::Vector{Float64}, Y::Matrix{Float64},
                    gopt::Vector{Float64})

Updates the gradient of the model in the case that the best optimal so far
'xopt' is not the origin of the sample set 'xbase'.

    - 'n': dimension of the search space.

    - 'm': number of interpolation conditions (points).

    - 'r': number of functions that make up the objective function fmin.

    - 'countf': function evaluation counter.

    - 'xopt': n-dimensional vector (optimal point so far).

    - 'hq': (n(n+1)/2)-dimensional vector (explicit second derivatives of the
    quadratic model).

    - 'pq': m-dimensional vector (parameters of the implicit second derivatives
    of the quadratic model).

    - 'Y': (n x m)-dimensional matrix (set of sample points).

The function modifies the argument:

    - 'gopt': n-dimensional vector (gradient of the quadratic model at 
    'xbase + xopt').

"""
function update_gopt!(
                        n::Int64,
                        m::Int64,
                        r::Int64,
                        countf::Int64,
                        xopt::Vector{Float64},
                        hq::Vector{Float64},
                        pq::Vector{Float64},
                        Y::Matrix{Float64},
                        gopt::Vector{Float64}
                        )

    idx_h = 0

    for j=1:n
        for i=1:j
            idx_h += 1
            if i < j
               gopt[j] += hq[idx_h] * xopt[i]
            end
            gopt[i] += hq[idx_h] * xopt[j]
        end
    end
    if countf > (r + m - 1)
        for k=1:m
            tmp = 0.0
            for j=1:n
                tmp += Y[j, k] * xopt[j]
            end
            tmp *= pq[k]
            for i=1:n
                gopt[i] += tmp * Y[i, k]
            end
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
    ρ = ( model.fval[model.kopt[]] - f_y ) / ( diff )

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
    ρ = ( model.fval[model.kopt[]] - fmin_y ) / ( diff )

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
                            fopt::Float64
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

        print_iteration( full_calc, it_flag, countit, countf, model.imin[], δ, Δ, model.fval[model.kopt[]] )

        if verbose ≥ 2

            print_warning(exit_flag)

            if model.kopt[] != 1

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
                                active_set::Vcetor{Bool}, 
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