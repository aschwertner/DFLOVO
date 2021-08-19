#-------------------------------------------------------------------------------
# Set of useful functions related to LinearModel struct.
#-------------------------------------------------------------------------------

function create_linear_model(n)

    return LinearModel(n, n + 1, Ref{Int64}(), Ref{Int64}(), Ref{Float64}(), zeros(Float64, n), zeros(Float64, n), zeros(Float64, n), zeros(Float64, n + 1), zeros(Float64, n), zeros(Float64, n, n))

end

function LinearModel(
                        n::Int64,
                        imin::Int64,
                        δ::Float64,
                        xbase::Vector{Float64},
                        fval::Vector{Float64},
                        Y::Matrix{Float64}
                        )

    dst = δ * ones(Float64, n)
    
    g = A \ ( fval[2:end] .- fval[1] )
    c = fval[1] - dot( g, xbase )

    return LinearModel(n, imin, c, g, xbase, fval, dst, Y)

end

( model::LinearModel )( x::AbstractVector ) = model.c + dot( model.g, x )

∇( model::LinearModel )( x::AbstractVector ) = model.g

∇2( model::LinearModel )( x::AbstractVector ) = false * I

function update_model!(
                        model::LinearModel, 
                        t::Int64,
                        fval_d::Float64,
                        d::Vector{Float64}
                        )

    # Updates Y set.
    for i = 1:model.n
        @. model.Y[i, :] -= d
    end

    # Puts the old xbase in the t-th vector in Y.
    @. model.Y[t, :] =  - d
    model.fval[t + 1] = model.fval[1]

    # Updates the distance vector.
    for i = 1:model.n
        model.dst[i] = norm(model.Y[i, :])
    end

    # Adds the new center.
    model.fval[1] = fval_d
    @. model.xbase += d 

    return rebuild_model!(model)

end

function model_from_scratch!(
                                model::LinearModel,
                                func_list::Array{Function, 1}, 
                                δ::Float64,
                                b::Vector{Float64}
                                )

    # Sets Y to a empty matrix.
    @. model.Y = 0.0

    for i = 1:model.n
        
        # Computes the new interpolation points and the distances.
        α = min( b[i] - model.xbase[i], δ )
        model.Y[i, i] = α
        model.dst[i] = α

        # Computes the function values.
        model.xbase[i] += α
        model.fval[i + 1] = func_list[model.imin](model.xbase)
        model.xbase[i] -= α

    end

    return rebuild_model!(model)

end

function rebuild_model!(
                        model::LinearModel
                        )

    # Computes the new c and g.
    # !!! Very inefficient !!!

    model.g = model.Y \ ( model.fval[2:end] .- model.fval[1] )
    model.c = model.fval[1] - dot( model.g, model.xbase )

    return LinearModel(model.n, model.imin, model.c, model.g, model.xbase, model.fval, model.dst, model.Y )

end

function active_set!(
                        model::LinearModel,
                        ao::Vector{Float64},
                        bo::Vector{Float64},
                        active_idx::Vector{Bool}
                        )
    for i=1:model.n
        if ( model.Y[model.kopt, i] == ao[i] && model.g[i] >= 0.0 ) || ( model.Y[model.kopt, i] == bo[i] && model.g[i] <= 0.0 )
            active_idx[i] = true
        else
            active_idx[i] = false
        end
    end

end

function trsbox!(
                    model::LinearModel,
                    Δ::Float64,
                    ao::Vector{Float64},
                    bo::Vector{Float64},
                    active_idx::Vector{Bool},
                    x::Vector{Float64},
                    d::Vector{Float64}
                    )
    
    # Copies the shifted vector 'xopt' to 'x'.
    copyto!(x, model.Y[model.kopt, :])

    # Creates the active constraint vector 'active_idx'.
    active_set!(model, ao, bo, active_idx)

    # Computes the symmetric vector of the projection of the 
    # gradient of the model by the set of active constraints.
    projection_active_set!(model.g, active_idx, d, sym = true)

    # If the set of active constraints is not complete,
    # calculates step 'd' and updates point 'x'.
    if sum(active_idx) != model.n
        α = Inf
        α_B = Inf
        α_Δ = Δ / norm(d)

        for i=1:model.n
            if d[i] > 0.0
                α = ( bo[i] - model.Y[model.kopt, i] ) / d[i]
            elseif d[i] < 0.0
                α = ( ao[i] - model.Y[model.kopt, i] ) / d[i]
            end
            if α < α_B
                α_B = α
            end
        end

        if α_Δ ≤ α_B
            @. d *= α_Δ
        else
            @. d *= α_B
        end

        @. x += d

    end
    
end