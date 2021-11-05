#-------------------------------------------------------------------------------
# Set of useful functions related to LinearModel struct.
#-------------------------------------------------------------------------------

export QuadraticModel

"""

Define the type for quadratic models.

    - 'n': model dimension.
    - 'm': number of interpolation points.
    - 'imin': index of the function interpolated.
    - 'kopt': index of the best point so far in 'Y'.
    - 'g': gradient of the model.
    - 'gopt': gradient of the model in 'xopt'.
    - 'hq': holds the explicit second derivatives of the quadratic model.
    - 'pq': holds parameters of the implicit second derivatives of the quadratic model.
    - 'xbase': origin of the sample set and center of the model.
    - 'xopt': best point so far (in terms of function values).
    - 'fval': set of the function values of the interpolation points.
    - 'dst': distances between 'xopt' and other interpolation points.
    - 'Y': set of interpolation points, shifted from the center of the sample set 'xbase'.
    - 'BMAT': holds the elements of 'Ξ', with the exception of its first column,
            and the elements of 'Υ', with the exception of its first row and column.
    - 'ZMAT': holds the elements of 'Z', from the factorization 'Ω = ZZ^T'.

"""
struct QuadraticModel <: AbstractModel

    n     :: Int64
    m     :: Int64
    imin  :: Ref{Int64}
    kopt  :: Ref{Int64}
    g     :: AbstractVector
    gopt  :: AbstractVector
    hq    :: AbstractVector
    pq    :: AbstractVector
    xbase :: AbstractVector
    xopt  :: AbstractVector
    fval  :: AbstractVector
    dst   :: AbstractVector
    Y     :: AbstractMatrix
    BMAT  :: AbstractMatrix
    ZMAT  :: AbstractMatrix

end


"""
    
    QuadraticModel(n::Int64, m::Int64)

Create a quadratic interpolation model for a function (not specified yet) of R^{n}, with 'm' interpolation points.

    - 'n': dimension of the model.
    - 'm': number of interpolation points.

Returns a model of QuadraticModel type.

"""
QuadraticModel(n::Int64, m::Int64) = QuadraticModel(
    n, m, Ref{Int64}(), Ref{Int64}(), 
    zeros(Float64, n), zeros(Float64, n), zeros(Float64, convert(Int64, n * ( n + 1 ) / 2) ), 
    zeros(Float64, m), zeros(Float64, n), zeros(Float64, n), zeros(Int64, m), 
    zeros(Float64, m - 1), zeros(Float64, m - 1, n), zeros(Float64, n + m, n), zeros(Float64, m, m - n - 1) )

"""

    construct_model!(func_list::Array{Function, 1}, imin_idx::Int64, 
                        δ::Float64, fbase::Float64, xbase::Vector{Float64},
                        ao::Vector{Float64}, bo::Vector{Float64},
                        model::QuadraticModel)

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

    - 'model': model of QuadraticModel type.

"""
function construct_model!(
                            func_list::Array{Function, 1},
                            imin_idx::Int64, 
                            δ::Float64, 
                            fbase::Float64,
                            xbase::Vector{Float64},
                            ao::Vector{Float64},
                            bo::Vector{Float64},
                            model::QuadraticModel
                            )

    # Saves information about the index i ∈ I_{min}('xbase') and f_{i}('xbase').
    model.imin[] = imin_idx
    model.fval[1] = fbase
    
    # Saves information about 'model.xbase'.
    @. model.xbase = xbase

    # Usefull constant and variables
    aux = 1.0 / δ ^ 2.0
    kopt = 0
    fopt = fbase

    for i = 1:( model.m - 1 )

        # Usefull variables
        j = i - model.n
        k = i + 1

        if i ≤ ( 2 * model.n )

            if i ≤ model.n

                if bo[i] == 0.0

                    α = - δ

                else

                    α = δ

                end

                model.Y[i, i] = α
                model.dst[i] = δ

            else

                α = model.Y[j, j]

                if ao[j] == 0.0

                    β = min( 2.0 * δ, bo[j] )
                
                elseif bo[j] == 0.0

                    β = max( - 2.0 * δ, ao[j] )
                
                else

                    β = - δ

                end

                model.Y[i, j] = β
                model.dst[i] = abs( β )

            end

        else

            tmp = div( i - model.n - 1, model.n )
            p_j = i - tmp * model.n - model.n
            q_j = p_j + tmp

            if q_j > model.n

                tmp = p_j
                p_j = q_j - model.n
                q_j = tmp

            end

            model.Y[i, p_j] = model.Y[p_j, p_j]
            model.Y[i, q_j] = model.Y[q_j, q_j]
            model.dst[i] = sqrt( model.Y[i, p_j] ^ 2.0 + model.Y[i, q_j] ^ 2.0 )

        end

        # Uses 'xbase' has workspace to reconstruct the sample points.
        reconstruct_original_point!( model, i, a, b, ao, bo, xbase )
        model.fval[k] = fi_eval( func_list, imin_idx, xbase )

        # Searches for the least function value to determine 'kopt'.
        if model.fval[k] < fopt

            kopt = i
            fopt = model.fval[k]

        end

        # Determines 'g', 'hq', 'BMAT' and 'ZMAT'. 
        if i ≤ ( 2 * model.n )

            if i ≤ model.n

                model.g[i] = ( model.fval[k] - fbase ) / α

                if model.m < ( k + model.n )

                    model.BMAT[1, i] = - 1.0 / α
                    model.BMAT[k, i] = 1.0 / α
                    model.BMAT[model.m + i, i] = - 0.5 * δ ^ 2.0

                end

            else

                idx_h = convert( Int64, j * ( j + 1 ) / 2 )
                tmp = ( model.fval[k] - fbase ) / β
                diff = β - α
                model.hq[idx_h] = 2.0 * ( tmp - model.g[j] ) / diff
                model.g[j] = ( model.g[j] * β - tmp * α ) / diff

                # Reorders the sample points if necessary.
                if α * β < 0.0

                    if model.fval[k] < model.fval[k - model.n]

                        tmp = model.fval[k - model.n]
                        model.fval[k - model.n] = model.fval[k]
                        model.fval[k] = tmp

                        if kopt == i

                            kopt = i - model.n

                        end

                        model.Y[j, j] = β
                        model.Y[i, j] = α
                        model.dst[j] = abs( β )
                        model.dst[i] = abs( α )

                    end

                end

                model.BMAT[1, j] = - ( α + β ) / ( α * β )
                model.BMAT[k, j] = - 0.5 / model.Y[i - model.n, j]
                model.BMAT[k - model.n, j] = - model.BMAT[1, j] - model.BMAT[k, j]
                model.ZMAT[1, j] = sqrt( 2.0 ) / ( α * β )
                model.ZMAT[k, j] = sqrt( 0.5 ) / δ ^ 2.0
                model.ZMAT[k - model.n, j] = - model.ZMAT[1, j] - model.ZMAT[k, j]

            end

        else

            idx_h = convert( Int64, q_j * ( q_j - 1 ) / 2 + p_j ) 
            model.ZMAT[1, j] = aux
            model.ZMAT[k, j] = aux
            model.ZMAT[q_j + 1, j] = - aux
            model.ZMAT[p_j + 1, j] = - aux
            tmp = model.Y[i, q_j] * model.Y[i, p_j]
            model.hq[idx_h] = ( fbase - model.fval[q_j + 1] - model.fval[p_j + 1] + model.fval[k] ) / tmp

        end

    end

    # Reconstructs the point 'xopt'
    reconstruct_original_point!( model, kopt, a, b, ao, bo, xbase )
    @. model.xopt = xbase
    model.kopt[] = kopt

    # Computes 'gopt' and updates 'dst' if necessary
    if kopt == 0

        @. model.gopt = model.g

    else

        update_gopt!( model, true )

        for i=1:(model.m - 1)

            if i == kopt

                continue

            else

                model.dst[i] = norm( model.Y[i, :] - model.Y[kopt, :] )

            end

        end

    end

end

"""

    update_gopt!(model::QuadraticModel, first_call::Bool)

Computes the gradient of the model in 'xopt'.

    - 'model': model of QuadraticModel type.
    
    - 'first_call' boolean that indicates if the model was newly created.

"""
function update_gopt!( 
                        model::QuadraticModel,
                        first_call::Bool
                        )
    
    # Computes H*d, where d = xopt - xbase, and store in 'gopt'.
    mul_hess_vec!( model, model.Y[model.kopt[], :], model.gopt )
    
    # Computes gopt = g + H * d
    @. model.gopt += model.g

    # Why? Search in BOBYQA's paper and implementation
    if !( first_call )

        for k=1:( model.m - 1 )

            tmp = model.pq[k] * dot( model.Y[k, :], model.Y[model.kopt[], :] )

            for i=1:model.n

                model.gopt[i] += tmp * Y[k, i]

            end

        end

    end

end

"""

    mul_hess_vec!(model::QuadraticModel, vec::Vector{Float64}, x::Vector{Float64})

Calculates the product between the model's Hessian and a vector, i.e., computes H * vec.

    - 'model': model of QuadraticModel type.

    - 'vec': n-dimensional vector.

The function modifies the argument:

    - 'x': n-dimensional vector (solution)

"""
function mul_hess_vec!(
                        model::QuadraticModel,
                        vec::Vector{Float64},
                        x::Vector{Float64}
                        )

    
    for j=1:model.n

        x[j] = 0.0

        for i=1:j

            idx_h += 1

            if i < j

                x[j] += model.hq[idx_h] * vec[i]

            end

            x[i] += model.hq[idx_h] * vec[j]

        end

    end

end

"""

    mul_vec_hess_vec!(model::QuadraticModel, v::Vector{Float64}, w::Vector{Float64})

Computes the product v^{T} * H * w, where H is the hessian of the model.

    - 'model': model of QuadraticModel type.

    - 'v': n-dimensional vector.

    - 'w': n-dimensional vector.

Return the value of the product.

"""
function mul_vec_hess_vec!(
                            model::QuadraticModel,
                            v::Vector{Float64},
                            w::Vector{Float64}
                            )

    aux = 0.0

    for j=1:model.n

        for i=1:j

            idx_h += 1

            if i < j

                aux += model.hq[idx_h] * w[i] * v[j]

            end

            aux += model.hq[idx_h] * w[j] * v[i]

        end

    end

    return aux

end

"""

    stationarity_measure(model::LinearModel, a::Vector{Float64}, 
                            b::Vector{Float64})

Calculates the stationarity measure π_{k} = || P_{Ω}( x_{k} - gopt_{k} ) - x_{k} ||.

    - 'model': model of QuadraticModel type. 

    - 'a': n-dimensional vector with the lower bounds.

    - 'b': n-dimensional vector with the upper bounds.

Returns the stationarity measure.

"""
function stationarity_measure(
                                model::QuadraticModel,
                                a::Vector{Float64}, 
                                b::Vector{Float64}
                                )

    aux = 0.0
    for i = 1:model.n

        aux += ( min( max( a[i], model.xopt[i] - model.gopt[i] ), b[i] ) - model.xopt[i] ) ^ 2.0

    end

    return sqrt(aux)

end

"""

    compute_active_set!(model::LinearModel, a::Vector{Float64},
                        b::Vector{Float64}, active_set::Vector{Bool})

Computes the set of active constraints.

    - 'model': model of QuadraticModel type.

    - 'a': n-dimensional vector with the lower bounds.

    - 'b': n-dimensional vector with the upper bounds.

The function modifies the argument:

    - 'active_set': boolean n-dimensional vector with the active constraints.

"""
function compute_active_set!(
                        model::QuadraticModel,
                        a::Vector{Float64},
                        b::Vector{Float64},
                        active_set::Vector{Bool}
                        )

    for i=1:model.n

        if ( model.xopt[i] == a[i] && model.gopt[i] ≥ 0.0 ) || ( model.xopt[i] == b[i] && model.gopt[i] ≤ 0.0 )

            active_set[i] = true

        else
            
            active_set[i] = false

        end

    end

end