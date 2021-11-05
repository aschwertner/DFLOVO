#-------------------------------------------------------------------------------
# Set of useful functions related to LinearModel struct.
#-------------------------------------------------------------------------------

export LinearModel, trsbox!

"""

Define the type for linear models.

    - 'n': model dimension.
    - 'm': number of interpolation points.
    - 'imin': index of the function interpolated.
    - 'kopt': index of the best point so far in 'Y'.
    - 'c': constant of the model.
    - 'g': gradient of the model.
    - 'xbase': origin of the sample set and center of the model.
    - 'xopt': best point so far (in terms of function values).
    - 'fval': set of the function values of the interpolation points.
    - 'dst': distances between 'xopt' and other interpolation points.
    - 'τY': auxiliary vector of QR factorization.
    - 'jpvtY': auxiliary vector with the permutations of QR factorization.
    - 'factorsY': auxiliary matrix of QR factorization.
    - 'Y': set of interpolation points, shifted from the center of the sample set 'xbase'.

"""
struct LinearModel <: AbstractModel

    n        :: Int64
    m        :: Int64
    imin     :: Ref{Int64}
    kopt     :: Ref{Int64}
    c        :: Ref{Float64}
    g        :: AbstractVector
    xbase    :: AbstractVector
    xopt     :: AbstractVector
    fval     :: AbstractVector
    dst      :: AbstractVector
    τY       :: AbstractVector
    jpvtY    :: AbstractVector
    factorsY :: AbstractMatrix
    Y        :: AbstractMatrix
    
end

"""
    
    LinearModel(n::Int64)

Create a linear interpolation model for a function (not specified yet) of R^{n}.

    - 'n': dimension of the model.

Returns a model of LinearModel type.

"""
LinearModel(n::Int64) = LinearModel(
    n, n + 1, Ref{Int64}(), Ref{Int64}(), Ref{Float64}(), 
    zeros(Float64, n), zeros(Float64, n), zeros(Float64, n), zeros(Float64, n + 1),
    zeros(Float64, n), zeros(Float64, n), zeros(Int64, n), zeros(Float64, n, n), zeros(Float64, n, n))


( model::LinearModel )( x::AbstractVector ) = model.c[] + dot( model.g, x )

"""

    rebuild_model!(model::LinearModel)

Computes the QR factorization of the matrix 'Y', the gradient and the constant of the model.

The function modifies the argument:

    - 'model': model of LinearModel type.

"""
function rebuild_model!(
                        model::LinearModel
                        )

    # Computes the QR-factorization of the matrix Y and stores the information in 'model.factorsY', 'model.τY' and 'model.jpvtY'.
    @. model.factorsY = model.Y
    qrY = qr!(model.factorsY, Val(true))
    @. model.τY = qrY.τ
    @. model.jpvtY = qrY.jpvt

    # Solves the system to determine the gradient of the linear model 'model.g'.
    @. model.g = model.fval[2:end] - model.fval[1]
    ldiv!(qrY, model.g)

    # Determines the constant of the linear model 'model.c'
    model.c[] = model.fval[1] - dot(model.g, model.xbase)

end

"""

    construct_model!(func_list::Array{Function, 1}, imin_idx::Int64, 
                        δ::Float64, fbase::Float64, xbase::Vector{Float64},
                        ao::Vector{Float64}, bo::Vector{Float64},
                        model::LinearModel)

Constructs the model based in the given information.

    - 'func_list': list containing the functions that determine the objective
    function fmin.

    - 'imin_idx': index belonging to the set I_{min}('xbase').

    - 'δ': radius of the sample set.

    - 'fbase': objective function value in 'xbase'.

    - 'xbase': n-dimensional vector (origin of the sample set).

    - 'ao': n-dimensional vector with the shifted lower bounds.

    - 'bo': n-dimensional vector with the shifted upper bounds.

The function modifies the argument:

    - 'model': model of LinearModel type.

"""
function construct_model!(
                            func_list::Array{Function, 1},
                            imin_idx::Int64, 
                            δ::Float64, 
                            fbase::Float64,
                            xbase::Vector{Float64},
                            ao::Vector{Float64},
                            bo::Vector{Float64},
                            model::LinearModel
                            )

    # Saves information about the index i ∈ I_{min}('xbase') and f_{i}('xbase').
    model.imin[] = imin_idx
    model.fval[1] = fbase

    kopt = 0

    for i = 1:model.n

        # Saves information about 'model.xbase' and 'model.xopt'.
        model.xbase[i] = xbase[i]
        model.xopt[i] = xbase[i]

        # Constructs the set interpolation points 'model.Y', shifted from 'xbase'.
        # Evaluates the function values and saves the information in 'model.fval'.
        if bo[i] == 0.0

            model.Y[i, i] = - δ
            xbase[i] -= δ
            model.fval[i + 1] = fi_eval(func_list, imin_idx, xbase)
            xbase[i] += δ

        else

            model.Y[i, i] = δ
            xbase[i] += δ
            model.fval[i + 1] = fi_eval(func_list, imin_idx, xbase)
            xbase[i] -= δ

        end

        # Saves the distance between the i-th point in 'model.Y' and 'model.xbase'.
        model.dst[i] = δ

        # Searches for the least function value to determine 'kopt'.
        if model.fval[i + 1] < fbase

            kopt = i
            fbase = model.fval[i + 1]

        end

    end

    # Saves the index of the best point in 'model.kopt'
    model.kopt[] = kopt

    # If the best point is not 'xbase', updates the vector 'xopt' and saves the information in 'model.xopt'.
    # Also updates the vector of 'model.dst', to store the distance between the interpolation points and 'xopt'.
    # Note that in the 'kopt' position of the 'model.dst' vector we keep the distance between 'xbase' and 'xopt' points. 
    if kopt != 0

        model.xopt[kopt] += model.Y[kopt, kopt]
        @. model.dst *= sqrt(2.0)
        model.dst[kopt] = δ

    end

    # Computes the constant and gradient of the linear model and the QR-factorization of the system.
    rebuild_model!( model )

end

"""

    update_model!(t::Int64, fnew::Float64, xnew::Vector{Float64},
                    model::LinearModel; trsbox_step::Bool=false)

Updates the model based in the given information.

    - 't': index of the point that must leave the sample set.

    - 'fnew: funtion value of the nw point.

    - 'xnew': n-dimensional vector (new sample point).

The function modifies the argument:

    - 'model': model of LinearModel type.

The optional argument is:

    - 'trsbox_step': Indicates if the step is of trust-region type. Set to 'false' by default.

"""
function update_model!(
                        t::Int64,
                        fnew::Float64,
                        xnew::Vector{Float64},
                        model::LinearModel;
                        trsbox_step::Bool=false
                        )

    if trsbox_step

        # If the step is of the TRSBOX type, the center of the sample set 'xbase' is never changed. Note that 't' is different from 'model.kopt[]'.
        # Since such a step is a descent direction for the model, 'fnew' is always less than the value of the function in 'xopt'.
        # The information about 'xopt' is also useful for the new model, so we maintain 'xopt' in the interpolation set 'Y', overwrite the point 't' by 'xnew',
        # and update 'xopt' to 'xnew'.

        model.kopt[] = t
        model.fval[t + 1] = fnew

        for i=1:model.n

            model.Y[t, i] = xnew[i] - model.xbase[i]
            model.xopt[i] = xnew[i]

        end

        @views for i=1:model.n

            if i == t

                model.dst[i] = norm( model.Y[i, :] )

            else

                model.dst[i] = norm( model.Y[i, :] - model.Y[t, :] )

            end

        end

    else

        # On the other hand, if the step is of type ALTMOV, the center of the sample set 'xbase' can be dropped or not from 'Y'.
        # Also, the best point so far 'xopt' can be changed or not, although the information about 'xopt' is always used for the new model.
        # Note that in the case where the points 'xbase' and 'xopt' coincide, the index 't' will never be 0. 

        if t == 0

            @assert model.kopt[] != 0 "In ALTMOV iterations, if the index 't' is set to 0, so the index 'kopt' cannot be 0."

            # The point 'xbase' is dropped.

            if fnew < model.fval[ model.kopt[] + 1 ]

                # The point 'xnew' is set to the new 'xbase' and 'xopt' is updated to 'xnew'.
                # The old information about the point 'xopt' is used in the new model.

                for i=1:model.n

                    @. model.Y[i, :] += model.xbase - xnew

                end

                for i=1:model.n

                    model.xbase[i] = xnew[i]
                    model.xopt[i] = xnew[i]
                    model.dst[i] = norm( model.Y[i, :] )

                end

                model.fval[1] = fnew
                model.kopt[] = 0

            else

                # The point 'xopt' is set to be the new 'xbase' and 'xnew' is incorporated as a point in the set 'Y'.

                for i=1:model.n

                    if i == model.kopt[]

                        @. model.Y[i, :] = xnew - model.xopt

                    else

                        @. model.Y[i, :] += model.xbase - model.xopt

                    end

                end

                for i=1:model.n

                    model.xbase[i] = model.xopt[i]
                    model.dst[i] = norm( model.Y[i, :] )

                end

                model.fval[1] = model.fval[ model.kopt[] + 1 ]
                model.fval[ model.kopt[] + 1 ] = fnew
                model.kopt[] = 0

            end

        else

            # The point 'xbase' is maintained.

            model.fval[ t + 1 ] = fnew
            @. model.Y[t, :] = xnew - model.xbase
            
            if fnew < model.fval[ model.kopt[] + 1 ]

                # The point 'xopt' is updated to 'xnew'.

                model.dst[t] = norm( model.Y[t, :] )
                model.kopt[] = t

                for i=1:model.n

                    model.xopt[i] = xnew[i]

                    if i != t

                        model.dst[i] = norm( model.Y[i, :] - model.Y[t, :] )

                    end

                end

            else

                model.dst[t] = norm( xnew - model.xopt )

            end

        end

    end

    # Computes the constant and gradient of the linear model and the QR-factorization of the system.
    rebuild_model!( model )

end

"""

    construct_new_model!(func_list::Array{Function, 1}, imin_idx::Int64, δ::Float64, 
                            fbase::Float64, xbase::Vector{Float64}, a::Vector{Float64},
                            b::Vector{Float64}, ao::Vector{Float64}, bo::Vector{Float64},
                            model::LinearModel)

Constructs the model based in the given information.

    - 'func_list': list containing the functions that determine the objective
    function fmin.

    - 'imin_idx': index belonging to the set I_{min}('xbase').

    - 'δ': radius of the sample set.

    - 'fbase': objective function value in 'xbase'.

    - 'xbase': n-dimensional vector (origin of the sample set).

    - 'a': n-dimensional vector with the lower bounds.

    - 'b': n-dimensional vector with the upper bounds.

The function modifies the argument:

    - 'ao': n-dimensional vector with the shifted lower bounds.

    - 'bo': n-dimensional vector with the shifted upper bounds.

    - 'model': model of LinearModel type.

"""
function construct_new_model!(
                            func_list::Array{Function, 1},
                            imin_idx::Int64, 
                            δ::Float64, 
                            fbase::Float64,
                            xbase::Vector{Float64},
                            a::Vector{Float64},
                            b::Vector{Float64},
                            ao::Vector{Float64},
                            bo::Vector{Float64},
                            model::LinearModel
                            )

    # Sets Y to a empty matrix.
    @. model.Y = 0.0

    # Saves information about the index i ∈ I_{min}('xbase') and f_{i}('xbase').
    model.imin[] = imin_idx
    model.fval[1] = fbase

    kopt = 0

    for i = 1:model.n

        # Saves information about 'model.xbase' and 'model.xopt'.
        model.xbase[i] = xbase[i]
        model.xopt[i] = xbase[i]

        # Computes the new relative bounds 'ao' and 'bo'.
        ao[i] = a[i] - xbase[i]
        bo[i] = b[i] - xbase[i]

        # Computes the new interpolation points 'model.Y', shifted from 'xbase'.
        # Evaluates the function values and saves the information in 'model.fval'.
        if bo[i] > δ

            model.Y[i, i] += δ
            xbase[i] += δ
            model.fval[i + 1] = fi_eval(func_list, imin_idx, xbase)
            xbase[i] -= δ

        else

            model.Y[i, i] -= δ
            xbase[i] -= δ
            model.fval[i + 1] = fi_eval(func_list, imin_idx, xbase)
            xbase[i] += δ

        end

        # Saves the distance between the i-th point in 'model.Y' and 'model.xbase'.
        model.dst[i] = δ

        # Searches for the least function value to determine 'kopt'.
        if model.fval[i + 1] < fbase

            kopt = i
            fbase = model.fval[kopt]

        end

    end

    # Saves the index of the best point in 'model.kopt'
    model.kopt[] = kopt

    # If the best point is not 'xbase', updates the vector 'xopt' and saves the information in 'model.xopt'.
    # Also updates the vector of 'model.dst', to store the distance between the interpolation points and 'xopt'.
    # Note that in the 'kopt' position of the 'model.dst' vector we keep the distance between 'xbase' and 'xopt' points. 
    if kopt != 0

        model.xopt[ kopt ] += model.Y[ kopt, kopt ]
        @. model.dst *= sqrt(2.0)
        model.dst[kopt] = δ

    end

    # Computes the constant and gradient of the linear model and the QR-factorization of the system.
    rebuild_model!( model )

end

"""

    stationarity_measure(model::LinearModel, a::Vector{Float64}, 
                            b::Vector{Float64})

Calculates the stationarity measure π_{k} = || P_{Ω}( x_{k} - g_{k} ) - x_{k} ||.

    - 'model': model of LinearModel type. 

    - 'a': n-dimensional vector with the lower bounds.

    - 'b': n-dimensional vector with the upper bounds.

Returns the stationarity measure.

"""
function stationarity_measure(
                                model::LinearModel,
                                a::Vector{Float64}, 
                                b::Vector{Float64}
                                )

    aux = 0.0
    for i = 1:model.n

        aux += ( min( max( a[i], model.xopt[i] - model.g[i] ), b[i] ) - model.xopt[i] ) ^ 2.0

    end

    return sqrt(aux)

end

"""

    compute_active_set!(model::LinearModel, a::Vector{Float64},
                        b::Vector{Float64}, active_set::Vector{Bool})

Computes the set of active constraints.

    - 'model': model of LinearModel type.

    - 'a': n-dimensional vector with the lower bounds.

    - 'b': n-dimensional vector with the upper bounds.

The function modifies the argument:

    - 'active_set': boolean n-dimensional vector with the active constraints.

"""
function compute_active_set!(
                        model::LinearModel,
                        a::Vector{Float64},
                        b::Vector{Float64},
                        active_set::Vector{Bool}
                        )

    for i=1:model.n

        if ( model.xopt[i] == a[i] && model.g[i] ≥ 0.0 ) || ( model.xopt[i] == b[i] && model.g[i] ≤ 0.0 )

            active_set[i] = true

        else
            
            active_set[i] = false

        end

    end

end

"""

    trsbox!(model::LinearModel, Δ::Float64, a::Vector{Float64}, b::Vector{Float64},
                active_set::Vector{Bool}, x::Vector{Float64}, d::Vector{Float64},
                s::Vector{Float64})

Computes a descent direction for the model based on TRSBOX routine of BOBYQA.

    - 'model': model of LinearModel type.

    - 'Δ': trust-region radius.

    - 'a': n-dimensional vector with the lower bounds.

    - 'b': n-dimensional vector with the upper bounds.

The function modifies the argument:

    - 'active_set': boolean n-dimensional vector with the active constraints.

    - 'x': n-dimensional vector (candidate for new sampling point 'x': 'xopt' + 'd').

    - 'd': n-dimensional vector (descent direction).

    - 's': n-dimensional vector (auxiliary vector for workspace).

Returns a flag indicating the stopping criteria.

"""
function trsbox!(
                    model::LinearModel,
                    Δ::Float64,
                    a::Vector{Float64},
                    b::Vector{Float64},
                    active_set::Vector{Bool},
                    x::Vector{Float64},
                    d::Vector{Float64},
                    s::Vector{Float64}
                    )

    # Copies vector 'xopt' to 'x' and sets vectors 'd' and 's' to zero.
    for i = 1:model.n

        x[i] = model.xopt[i]
        d[i] = 0.0
        s[i] = 0.0

    end

    # Saves the best function value so far.
    fopt = model.fval[ model.kopt[] + 1 ]

    # Creates the active constraint vector 'active_set'.
    compute_active_set!(model, a, b, active_set)

    # If the set of active constraints is full, 'd' is returned as the null vector.
    if sum(active_set) == model.n

        return :full_active_set
    
    end

    # Computes the symmetric vector of the projection of the gradient of the model 
    # by the set of active constraints and stores in 's'.
    projection_active_set!(model.g, active_set, s, sym = true)

    # If 's' is a null vector, we won't get any decreases in that direction.
    # Therefore, 'd' is returned as the null vector.
    if norm(s) ≈ 0.0

        return :null_search_direction

    end

    while true

        # Computes the step 'α' along the direction 's'.
        α, chose_αΔ, idx_B = compute_alpha_linear(Δ, a, b, x, d, s)

        # Updates the search direction 'd'.
        @. d += α * s

        if chose_αΔ

            # If the trust-region boundary has been reached (case 'α' = 'α_Δ'),
            # prepares for an auxiliary procedure that tries to increase the reduction
            # of the obtained model value.

            while true

                # Uses 'x' as work vector to store projection of direction 'd' by a set of active constraints 'active_set'.
                projection_active_set!(d, active_set, x)

                # Computes P_{I}(∇model(xopt+d)) = P_{I}('model.g') and stores in 's'.
                projection_active_set!(model.g, active_set, s)

                pdTpd = dot( x, x )
                pdTpg = dot( x, s )
                pgTpg = dot( s, s )

                # Tests the stopping criteria.
                if ( pdTpd * pgTpg - pdTpg ^ 2.0 ) ≤ 1.0e-4 * ( fopt - model.c[] - dot( model.g, model.xopt ) - dot( model.g, d ) )

                    step_projection!(model, a, b, x, d)

                    return :trust_region_boundary

                else

                    # Otherwise, computes a new search direction 's', spanned by P_{I}(d) and P_{I}(∇model(xopt+d)),
                    # that satisfies ||s|| = ||P_{I}(d)||, s'P_{I}(d) = 0, and s'P_{I}(∇model(xopt+d)) < 0.0.
                    # Note that 'x' stores P_{I}(d) and 's' stores P_{I}(∇model(xopt+d)).
                    # After the call of 'new_search_direction', the desired direction is allocated in 's'.
                    new_search_direction!(pdTpd, pdTpg, pgTpg, x, s)

                    # Verifies if the direction 's' is null.
                    if norm(s) ≈ 0.0

                        step_projection!(model, a, b, x, d)
            
                        return :null_new_search_direction
                
                    end

                    # Computes the parameter θ that satisfies a ≤ xopt + d(θ) ≤ b, and g'd(θ) < 0.0.
                    # Note that the vector 's' holds d(θ) - d.
                    chose_θQ = compute_theta_linear!(model, a, b, d, x, s)

                    if norm(s) ≈ 0.0

                        step_projection!(model, a, b, x, d)

                        return :null_angle_search_direction

                    end

                    # Updates the direction 'd' to be d(θ).
                    @. d += s

                    # Updates 'index_set' with indices for which xopt + d is restricted by one of the bounds on the variables.
                    update_active_set!( a, b, model.xopt, d, active_set )
                    
                    if chose_θQ

                        # Computes P_{I}(∇model(xopt+d)) = P_{I}('model.g') and stores in 's'.
                        projection_active_set!(model.g, active_set, s)

                        # Verifies the stopping criteria.
                        if ( norm(s) * Δ ) ≤ 1.0e-2 * ( fopt - model.c[] - dot( model.g, model.xopt ) - dot( model.g, d ) )

                            step_projection!(model, a, b, x, d)

                            return :tired_trust_region

                        end

                    end

                end

            end

        else

            # Otherwise, the current line search is restricted by a bound constraint. If the stooping criteria is true, the computational effort
            # to perform a new iteration is not promising and the procedure is terminated.

            # Computes - P_{I}(∇model(xopt+d)) = - P_{I}('model.g') and stores in 's'.
            projection_active_set!(model.g, active_set, s, sym = true)

            # Verifies the stopping criteria.
            if ( norm(s) * Δ ) ≤ 1.0e-2 * ( fopt - model.c[] - dot( model.g, x ) - dot( model.g, d ) )

                step_projection!(model, a, b, x, d)
                
                return :tired_bound_constraint

            else

                # Adds the index of the achieved bound constraint to 'active_set',
                # and updates 's' to be the symmetric projection of 'model.g' in the new set.

                if idx_B != 0

                    active_set[ idx_B ] = true
                    s[ idx_B ] = 0.0

                end
                
            end

        end

        # If the set of active constraints is full, the procedure terminates.
        if sum(active_set) == model.n

            step_projection!(model, a, b, x, d)

            return :full_active_set
    
        end

        # If the direction 's' is null, the procedure also terminates.
        if norm(s) ≈ 0.0

            step_projection!(model, a, b, x, d)

            return :null_search_direction
        
        end

    end
    
end

"""

    choose_index_trsbox(model::LinearModel, xnew::Vector{Float64}, v::Vector{Float64})

Chooses the point that must leave the sampling set for trust-region steps.

    - 'model': model of LinearModel type.

    - 'xnew': n-dimensional vector (new sample point).

The function modifies the argument:
    
    - 's': n-dimensional vector (auxiliary vector for workspace).

Returns an index 't' which indicates the point that must leave the sample set.

"""
function choose_index_trsbox(
                                model::LinearModel,
                                xnew::Vector{Float64},
                                v::Vector{Float64}
                                )
    
    qrY = QRPivoted( model.factorsY, model.τY, model.jpvtY )
    α = - Inf
    t = - 1

    for i = 1:model.n

        @. v = 0.0
        v[i] = 1.0
        ldiv!(qrY, v)
        α_t = abs( 1.0 + dot( xnew, v ) - dot( model.xbase, v ) - dot( model.Y[i, :], v ) )

        if α_t > α

            α = α_t
            t = i

        end

    end

    return t

end


"""

    lagrange_pol_t(model::LinearModel, idx_t::Int64, ∇Λ_t::Vector{Float64}, y::Vector{Float64})

Computes the t-th Lagrange polynomial at 'y', i.e., computes Λ_{t}(y).

    - 'model': model of LinearModel type.

    - 'idx_t': index of the Lagrange polynomial.

    - '∇Λ_t': n-dimensional vector (gradient of the t-th Lagrange polynomial).

    - 'y': n-dimensional vector (point of interest).

Returns the function value Λ_{t}(y).

"""
function lagrange_pol_t(
                        model::LinearModel,
                        idx_t::Int64,
                        ∇Λ_t::Vector{Float64},
                        y::Vector{Float64}
                        )

    if idx_t == 0

        return 1.0 + dot( y, ∇Λ_t ) - dot( model.xbase, ∇Λ_t )

    else

        return dot( y, ∇Λ_t ) - dot( model.xbase, ∇Λ_t )

    end

end

"""

    altmov!(model::LinearModel, idx_t::Int64, Δ::Float64, a::Vector{Float64}, 
                b::Vector{Float64}, x::Vector{Float64}, d::Vector{Float64},
                v::Vector{Float64}, w::Vector{Float64}, altmov_set::Vector{Bool})

Calculates a direction to improve the poisedness of the sample set based on ALTMOV routine of BOBYQA.

    - 'model': model of LinearModel type.

    - 'idx_t': index of the point that must leave the sample set.

    - 'Δ': trust-region radius.

    - 'a': n-dimensional vector with the lower bounds.

    - 'b': n-dimensional vector with the upper bounds.

The function modifies the argument:

    - 'x': n-dimensional vector (candidate for new sampling point 'x': 'xopt' + 'd').

    - 'd': n-dimensional vector (descent direction).

    - 'v': n-dimensional vector (auxiliary vector for workspace).

    - 'w': n-dimensional vector (auxiliary vector for workspace).

    - 'altmov_set': boolean n-dimensional vector with the active constraints.

Returns a flag indicating the stopping criteria.

"""
function altmov!(
                    model::LinearModel,
                    idx_t::Int64,
                    Δ::Float64,
                    a::Vector{Float64},
                    b::Vector{Float64},
                    x::Vector{Float64},
                    d::Vector{Float64},
                    v::Vector{Float64},
                    w::Vector{Float64},
                    altmov_set::Vector{Bool}
                    )

    # Reconstruct the QR-factorization of the system.
    qrY = QRPivoted( model.factorsY, model.τY, model.jpvtY )
    
    # Initializes some usefull variables and constants
    best_abs_ϕ = - 1.0
    best_idx = 0
    best_α = Inf
    a_0 = - Δ ^ 2.0
    a_1 = 0.0

    # Computes the gradient of the "idx_t"-th Lagrange polynomial and stores the information in vector 'v'.
    if idx_t == 0.0

        @. x = - 1.0
        ldiv!( v, qrY, x )

    else

        @. x = 0.0
        x[idx_t] = 1.0
        ldiv!( v, qrY, x )

    end

    # Calculates the "Usual" Alternative Step
    for j = 1:model.n

        # Uses 'd' as a workspace to save the difference between the j-th interpolation point and 'xopt'.
        # If 'xopt' and 'xbase' coincide, such differences are already calculated and available in the matrix 'model.Y'.
        # Otherwise, if 'xopt' is not the center of the sample set, and 'j' = 'kopt', 'd' holds the difference between 'xbase' and 'xopt'.
        if model.kopt[] == 0

            @. d = model.Y[ j, : ]

        else

            if j == model.kopt[]

                @. d = - model.Y[ model.kopt[], : ]

            else

                @. d = model.Y[ j, : ] - model.Y[ model.kopt[], : ]

            end

        end

        # Gets the bounds for α_{j} relative to the trust-region.
        roots = []
        a_2 = dot( d, d )
        solve_quadratic!(a_2, a_1, a_0, roots)
        
        α_upper = maximum(roots)

        if ( α_upper == NaN )

            for i = 1:model.n

                x[i] = model.xopt[i]
                d[i] = 0.0

            end

            return :invalid_steplenght

        else

            α_lower = minimum(roots)

        end

        # Adjusts the bounds to the box constraints.
        for i = 1:model.n

            if d[i] > 0.0

                α_lower = max( α_lower, ( a[i] - model.xopt[i] ) / d[i] )
                α_upper = min( α_upper, ( b[i] - model.xopt[i] ) / d[i] )

            elseif d[i] < 0.0

                α_lower = max( α_lower, ( b[i] - model.xopt[i] ) / d[i] )
                α_upper = min( α_upper, ( a[i] - model.xopt[i] ) / d[i] )

            end

        end

        # Since the Lagrange polynomial is linear, the optimal α is either lower or upper value.
        # Uses 'x' as a workspace to hold 'xopt + α * (y_{j} - xopt)'
        @. x = model.xopt + α_lower * d
        ϕ_lower = lagrange_pol_t( model, idx_t, v, x )

        @. x = model.xopt + α_upper * d
        ϕ_upper = lagrange_pol_t( model, idx_t, v, x )

        if abs( ϕ_lower ) > abs( ϕ_upper )

            α_j = α_lower
            abs_ϕ_j = abs( ϕ_lower )

        else

            α_j = α_upper
            abs_ϕ_j = abs( ϕ_upper )

        end

        # Compares the best ϕ_j value so far
        if abs_ϕ_j > best_abs_ϕ

            best_abs_ϕ = abs_ϕ_j
            best_idx = j
            best_α = α_j

        end

    end

    # Calculates the "Cauchy" Alternative Step

    Λ_1 = altmov_cauchy!(model, idx_t, Δ, a, b, v, d, w, altmov_set)
    Λ_2 = altmov_cauchy!(model, idx_t, Δ, a, b, v, x, w, altmov_set, sym = true)

    # Chooses the best step
    ( val, idx ) = findmax( [ best_abs_ϕ, abs( Λ_1 ), abs( Λ_2 ) ] )

    if idx == 1

        if model.kopt[] == 0

            @. d = best_α * model.Y[ best_idx, : ]
        
        else

            if best_idx == model.kopt[]

                @. d = - best_α * model.Y[ best_idx, :]

            else

                @. d = best_α * ( model.Y[ best_idx, : ] - model.Y[ model.kopt[], : ] )

            end

        end
        
        step_projection!(model, a, b, x, d)

        return :usual_altmov

    elseif idx == 2

        step_projection!(model, a, b, x, d)

    else

        @. d = x
        step_projection!(model, a, b, x, d)

    end

    return :cauchy_altmov

end

"""

    altmov_cauchy!(model::LinearModel, idx_t::Int64, Δ::Float64, a::Vector{Float64}, 
                    b::Vector{Float64}, ∇Λ_t::Vector{Float64}, s::Vector{Float64},
                    z::Vector{Float64}, active_set::Vector{Bool}; sym::Bool=false)

Calculates a direction to improve the poisedness of the sample set based on the Cauchy step of ALTMOV routine from BOBYQA.

    - 'model': model of LinearModel type.

    - 'idx_t': index of the point that must leave the sample set.

    - 'Δ': trust-region radius.

    - 'a': n-dimensional vector with the lower bounds.

    - 'b': n-dimensional vector with the upper bounds.

    - '∇Λ_t': n-dimensional vector (gradient of the t-th Lagrange polynomial).

The function modifies the argument:

    - 's': n-dimensional vector (direction).

    - 'z': n-dimensional vector (candidate for new sampling point 'x': 'xopt' + 's').

    - 'active_set': boolean n-dimensional vector with the active constraints.

The optional argument is:

    - 'sym': indicates if it is calculated based on the symetric direction of '∇Λ_t'. Set to 'false' by default.

Returns the function value of the t-th Lagrange polynomial in 'z', i.e., Λ_t(z).

"""
function altmov_cauchy!(
                        model::LinearModel,
                        idx_t::Int64,
                        Δ::Float64,
                        a::Vector{Float64},
                        b::Vector{Float64},
                        ∇Λ_t::Vector{Float64},
                        s::Vector{Float64},
                        z::Vector{Float64},
                        active_set::Vector{Bool};
                        sym::Bool=false
                        )

    if sym

        p = 1.0

    else

        p = 2.0

    end

    norm2_s = 0.0

    for i = 1:model.n

        if ( (- 1.0) ^ p ) * ∇Λ_t[i] > 0.0

            s[i] = a[i] - model.xopt[i]
            norm2_s += s[i] ^ 2.0

        elseif ( (- 1.0) ^ p ) * ∇Λ_t[i] < 0.0

            s[i] = b[i] - model.xopt[i]
            norm2_s += s[i] ^ 2.0

        else

            s[i] == 0.0

        end

    end

    if norm2_s ≤ ( Δ ^ 2.0 )

        @. z = model.xopt + s

        return lagrange_pol_t( model, idx_t, ∇Λ_t, z )

    end

    # Saves a copy of 's' in 'z'.
    copyto!( z, s )
    
    sum_free = 0.0
    sum_fixed = 0.0

    for i = 1:model.n

        if s[i] ≈ 0.0

            active_set[i] = true

        else

            sum_free += ∇Λ_t[i] ^ 2.0
            active_set[i] = false

        end

    end

    if sum_free ≈ 0.0

        s .= 0.0
        return 0.0

    end

    while true

        num_μ = Δ ^ 2.0 - sum_fixed

        if num_μ > 0.0

            μ = sqrt( num_μ / sum_free )

            stop = true
            sum_free = 0.0
            sum_fixed = 0.0

            for i = 1:model.n

                if active_set[i]

                    s[i] = z[i]
                    sum_fixed += z[i] ^ 2.0

                else

                    s[i] = ( (- 1.0) ^ ( p + 1 ) ) * μ * ∇Λ_t[i]

                    step_i = s[i] + model.xopt[i]

                    if ( step_i < a[i] ) || ( step_i > b[i] )
                    
                        stop = false
                        active_set[i] = true
                        sum_fixed += z[i] ^ 2.0

                    else

                        sum_free += ∇Λ_t[i] ^ 2.0

                    end

                end
                
            end

            if stop || ( sum_free ≈ 0.0 )

                break

            end

        else

            break

        end

    end

    @. z = model.xopt + s

    return lagrange_pol_t( model, idx_t, ∇Λ_t, z )

end

"""
    compute_theta_linear!(model::LinearModel, a::Vector{Float64}, b::Vector{Float64},
                            d::Vector{Float64}, proj_d::Vector{Float64}, s::Vector{Float64})

Calculates the parameter θ that satisfies 'a' ≤ 'xopt' + 'd(θ)' ≤ 'b', and 'g'd(θ)' < 0.0, in 'trsbox!' routine.

    - 'model': model of LinearModel type.

    - 'a': n-dimensional vector with the lower bounds.

    - 'b': n-dimensional vector with the upper bounds.

    - 'd': n-dimensional vector (direction).

    - 'proj_d': n-dimensional vector (projection of the direction 'd' on the set of active constraints).

The function modifies the argument:

    - 's': n-dimensional vector (new direction).

Returns 'true' if the parameter 'θQ' is chosen, and 'false' otherwise.

"""
function compute_theta_linear!(
                                model::LinearModel, 
                                a::Vector{Float64},
                                b::Vector{Float64},
                                d::Vector{Float64},
                                proj_d::Vector{Float64},
                                s::Vector{Float64}
                                )

    cond_B(θ) = cond_θB(θ, a, b, model.xopt, d, proj_d, s)

    θB = binary_search(0.0, 0.25 * π , cond_B, 1.0e-2)

    if θB == NaN

        @. s = 0.0

    else

        gTd = dot(model.g, model.g)
        gTpd = dot(model.g, proj_d)
        gTs = dot(model.g, s)

        cond_Q(θ) = cond_θQ_linear(θ, gTd, gTpd, gTs)

        θQ = binary_search(0.0, 0.25 * π , cond_Q, 1.0e-2)

        if θQ == NaN

            @. s = 0.0

        else

            θ = min(θB, θQ)

            @. s = - proj_d + cos(θ) * proj_d + sin(θ) * s

            if θ == θQ

                return true

            else

                return false

            end

        end

    end

end