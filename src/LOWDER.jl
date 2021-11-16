# Turn off precompilation during development
__precompile__(false)

#-------------------------------------------------------------------------------
# Main file implementing LOWDER in Julia
#-------------------------------------------------------------------------------

module LOWDER

    # Load dependencies

    using LinearAlgebra
    using Printf

    # Load code
    include("lowder_types.jl")
    include("linear_models.jl")
    include("lovo_utils.jl")
    include("utils.jl")
    include("lowder_output.jl")

    function lowder(
                    func_list::Array{Function,1},
                    x::Vector{Float64},
                    a::Vector{Float64}, 
                    b::Vector{Float64},
                    δ::Float64,
                    Δ::Float64;
                    m::Int64=(2 * length(x) + 1),
                    maxit::Int64=(1000 * length(x)),
                    maxfun::Int64=(1000 * length(func_list) * length(x)),
                    maxcrit::Int64=(m - 1),
                    nρmax::Int64=3,
                    Γmax::Int64=1,
                    verbose::Int64=0,
                    δmin::Float64=1.0e-4,
                    πmin::Float64=1.0e-4,
                    β::Float64=1.0,
                    τ1::Float64=0.6,
                    τ2::Float64=1.5,
                    τ3::Float64=2.0,
                    η::Float64=0.1,
                    η1::Float64=0.3,
                    η2::Float64=0.6,
                    filename::String=""
                    )

        # Calculates the search space dimension.
        n = length(x)

        # Calculates the number of functions that make up the objective function fmin
        r = length(func_list)

        # Verify algorithm initialization conditions.
        if m != ( n + 1 )
            m_min = n + 2
            m_max = convert(Int64, ( n + 1 ) * ( n + 2 ) / 2)
            @assert m_min ≤ m ≤ m_max "The number of interpolation points 'm' must satisfy m ∈ [$(m_min), $(m_max)] for quadratic models. In the case of linear models, 'm' must be defined as $(n + 1)."
        end
        @assert length(a) == n "The vector 'a' must have dimension $(n)."
        @assert length(b) == n "The vector 'b' must have dimension $(n)."
        @assert 0.0 < δ ≤ Δ "The radius of the sample set 'δ' must be positive and less or equal to the trust-region radius 'Δ'."
        @assert verify_initial_room(n, δ, a, b) "The radius of the initial sample set 'δ' is not suitable, it must satisfy a[i]- b[i] >= 2δ, for i = 1, ..., $(n)."
        @assert maxit > 0 "The parameter 'maxit' must be positive."
        @assert maxfun ≥ r "The parameter 'maxfun' must be greater than ou equal to $(r)."
        @assert Γmax > 0 "The parameter 'Γmax' must be positive."
        @assert δmin > 0.0 "The parameter 'δmin' must be positive."
        @assert β > 0.0 "The parameter 'β' must be positive."
        @assert 0.0 < τ1 < 1.0 "The parameter 'τ1' must be positive and less than one."
        @assert 1.0 ≤ τ2 "The parameter 'τ2' must be greater than or equal to one."
        @assert τ2 ≤ τ3 "The parameter 'τ3' must be greater than or equal to 'τ2'."
        @assert 0.0 < η1 < 1.0 "The parameter 'η1' must be positive and less than one."
        @assert 0 ≤ η < η1 "The parameter 'η' must be nonnegative and less than 'η1'."
        @assert η1 ≤ η2 "The parameter 'η2' must be greater than or equal to 'η1'."
        @assert 0 ≤ verbose ≤ 3 "The parameter 'verbose' must be 0, 1, 2 or 3."

        # Verifies if the usuary wants to save a file with information of the iterations.
        if filename == ""

            saveinfo = false

        else

            saveinfo = true

        end
        
        # Sets some useful constants.
        Δinit = Δ
        δcrit = exp( log( ( 1.0 - sqrt( eps( Float64 ) ) ) * δmin ) - maxcrit * log( τ1 ) )

        # Initializes useful variables, vectors, and matrices.
        nρ = 0                        # Auxiliary counter for simplified 'ρ' calculations.
        nΓ = 0                        # Auxiliary counter for radii adjustments phase.
        nit = 0                       # Counts the number of iterations.
        nf = 0                        # Counts the number of 'f_{i}' function evaluations.
        ncrit = 0                     # Counts the number of consecutive criticality iterations.
        naltmov = 0                   # Counts the number of consecutive iterations of type ALTMOV, ignoring criticality iterations.
        fi_x = NaN                    # New function value.
        pred_red = NaN                # Predicted reduction.
        real_red = NaN                # Real reduction.
        ρ = NaN                       # Relative reduction.
        π = NaN                       # Stationarity measure.
        norm_d = NaN                  # Norm of the computed direction 'd'.
        δold = NaN                    # Radius of the sample set at the beginning of the iteration.
        Δold = NaN                    # Radius of the trust-region at the beginning of the iteration.
        full_calc = false             # Indicates that ρ was calculated using fmin.
        verif_new_info = true         # Indicates the need to verify the information obtained in the last iteration. 
        ao = zeros(Float64, n)        # Difference between the lower bounds 'a' and the center of the sample set, given by 'xbase'.
        bo = zeros(Float64, n)        # Difference between the upper bounds 'b' and the center of the sample set, given by 'xbase'.
        d = zeros(Float64, n)         # TRSBOX or ALTMOV direction.
        aux_v = zeros(Float64, n)     # Auxiliar vector for workspace.
        aux_w = zeros(Float64, n)     # Auxiliar vector for workspace.
        active_set = zeros(Bool, n)   # Set of active constraints.
        imin_set = zeros(Bool, r)     # Set I_{min}(x)
        it_flag = :nonspecified       # Initializes the iteration flag.
        exit_flag = :nonspecified     # Initializes the exit flag.
        status_flag = :nonspecified   # Initializes the status flag.

        #-------------------- Preparations for the first iteration ---------------------

        # Create an empty structure called 'model', which can hold a linear or quadratic model.
        if n == (m - 1)

            model = LinearModel(n)

        end

        # Modifies the initial estimate 'x' to be suitable for building the first model. 
        # Modifies 'ao' and 'bo' to store the 'a-x' and 'b-x' differences, respectively.
        correct_guess_bounds!(n, δ, a, b, x, ao, bo)

        # Computes the value of f_{min}('x') and an index 'imin_idx' ∈ I_{min}('x').
        fi_x, imin_idx = fmin_eval!(func_list, r, x, imin_set)

        # Updates de function call counter.
        nf += r

        # Creates the initial model by modifying the structure of the previous 'model' 
        # and calculates the QR factorization of the matrix M associated with the model definition.
        construct_model!(func_list, imin_idx, δ, fi_x, x, ao, bo, model)

        # Updates de function call counter.
        nf += n - 1

        if saveinfo

            data_file = open(filename, "w")

        end

        # Returns if 'nf' exceeds 'maxfun'.
        if nf ≥ maxfun

            it_flag = :nonspecified
            exit_flag = :max_evaluations

            # Creates the LOWDEROutput.
            output = create_output(model, exit_flag, full_calc, nit, nf, δold, Δold, π)

            # Prints information about the iteration, exit flag and LOWDEROutput.
            print_info(model, output, it_flag, status_flag, full_calc, verbose, δ, Δ, π, ρ, pred_red, real_red, d)

            if saveinfo
    
                save_info!(model, it_flag, status_flag, full_calc, nit, nf, δ, Δ, π, ρ, pred_red, real_red, d, data_file)
                println(data_file, exit_flag)
                close(data_file)
            
            end

            return output
            
        end
        
        # Main loop
        while true

            it_flag = :nonspecified
            status_flag = :nonspecified

            # Saves old information about radii
            δold = δ
            Δold = Δ

            # Computes the stationarity measure π = ||P_{Ω}(xopt - g) - xopt||
            π = stationarity_measure(model, a, b)

            # Verifies if 'δ' and 'π' are less than or equal to 'δmin' and 'πmin', respectively.
            if ( ( δ ≤ δmin ) && ( π ≤ πmin ) )

                # Sets iteration and exit flags
                exit_flag = :success
                break

            end

            if δ > β * π

                #------------------------------ Criticality phase ------------------------------

                # Sets iteration flag
                it_flag = :criticality

                # Update parameters and counter
                ncrit += 1
                naltmov = 0
                ρ = NaN
                pred_red = NaN
                real_red = NaN
                norm_d = NaN

                if ncrit ≥ maxcrit

                    # It allows performing at most 'maxcrit' criticality iterations until the end of the algorithm. 
                    δ = min( δ, δcrit )
                    ncrit = 0

                else
                    
                    δ *= τ1

                end

            else
                
                # Update counter
                ncrit = 0

                #------------------------------- Step calculation ------------------------------

                status_flag = trsbox!(model, Δ, a, b, active_set, x, d, aux_v)
                norm_d = norm( d )

                # Verifies if 'trsbox' returns a unexpected direction.
                if norm_d == NaN

                    # Sets iteration and exit flags
                    it_flag = :trust_region
                    exit_flag = :numerical_error
                    break

                elseif norm_d ≥ 0.5 * Δ

                    #------------------------------- Step acceptance -------------------------------

                    mnew = model(x)
                    pred_red = model.fval[ model.kopt[] + 1 ] - mnew
                    full_calc = false

                    if pred_red < 0

                        it_flag = :nonspecified
                        exit_flag = :nondescent
                        break

                    else

                        if nρ < nρmax

                            ρ, real_red, fi_x, imin_idx = relative_reduction(model, func_list, pred_red, x)
                            nf +=1

                        else

                            full_calc = true
                            ρ, real_red, fi_x, imin_idx = relative_reduction!(model, func_list, r, pred_red, x, imin_set)
                            nf += r - 1
                            nρ = 0

                        end

                    end

                else

                    # Update parameters
                    ρ = 0.0
                    pred_red = 0.0
                    real_red = 0.0

                end

                #------------------------------- Radii updates ---------------------------------

                if ρ < η1

                    if ρ > 0.0

                        # Iterations of type 'trust_region' and 'bad_trust_region'.

                        Δ = max( τ1 * Δ ,  δ )

                    else

                        # Iterations of type 'altmov'.

                        if naltmov > ( model.m - 1 )

                            δ *= τ1
                            Δ *= τ1

                        else

                            Δ = max( τ1 * Δ ,  δ )

                        end

                    end

                elseif ( ρ > η2 ) && ( norm_d ≈ Δ )

                    δ *= τ2
                    Δ *= τ2

                end

            end

            # Chooses the point 't' that must leave the interpolation set. In the case of an ALTMOV call,
            # it also computes the new point 'x' and the direction 'd'.
            if ρ ≥ η

                # Sets iteration flag
                it_flag = :trust_region

                # Chooses the point that must leave the interpolation set 'model.Y'.
                t = choose_index_trsbox(model, x, aux_v)

                # Renitializes the counter
                nρ = 0
                naltmov = 0

            elseif ρ > 0.0

                # Sets iteration flag
                it_flag = :bad_trust_region

                # Chooses the point that must leave the interpolation set 'model.Y'.
                t = choose_index_trsbox(model, x, aux_v)

                # Updates the counters.
                nρ += 1
                naltmov = 0

            else

                # Sets iteration flag
                if it_flag != :criticality

                    it_flag = :altmov

                    if naltmov > ( model.m - 1 )

                        naltmov = 0

                    else

                        naltmov += 1

                    end

                end

                # Chooses the point that must leave the interpolation set 'model.Y'.
                t = choose_index_altmov(model)

                # Computes the new interpolation point.
                status_flag = altmov!(model, t, δold, a, b, x, d, aux_v, aux_w, active_set)

                # Verifies if 'altmov' returns a unexpected direction.
                if norm( d ) == NaN

                    exit_flag = :numerical_error
                    break

                end

                # Computes the new function value.
                fi_x = fi_eval(func_list, model.imin[], x)

                # Updates the counters.
                nf += 1

            end

            #------------------------- Verifies output conditions --------------------------
        
            # Verifies if 'nit' exceeds 'maxit', or if 'nf' exceeds 'maxfun'. 
            if nit ≥ maxit

                exit_flag = :max_iterations
                break

            elseif nf ≥ maxfun

                exit_flag = :max_evaluations
                break

            elseif δ < eps(Float64)

                exit_flag = :small_sampling_radius
                break

            else

                if verbose ≥ 1

                    print_iteration(it_flag, status_flag, full_calc, nit, nf, model.imin[], δold, Δold, π, ρ, pred_red, real_red, model.fval[ model.kopt[] + 1 ], model.xopt, d)

                end

                if saveinfo
    
                    save_info!(model, it_flag, status_flag, full_calc, nit, nf, δold, Δold, π, ρ, pred_red, real_red, d, data_file)
                
                end

            end

            #--------------------------- Radii adjustments phase ---------------------------

            if ρ > 0.0

                if full_calc

                    if imin_set[ model.imin[] ]

                        imin_idx = model.imin[]

                    end

                else

                    fi_x, imin_idx = fmin_partial_eval!( func_list, r, model.imin[], fi_x, x, imin_set )
                    nf += r - 1

                    if imin_set[ model.imin[] ]

                        imin_idx = model.imin[]

                    end

                end

                if imin_idx != model.imin[]

                    if ρ ≥ η1

                        nΓ = 0

                    end

                    if nΓ ≤ Γmax

                        # Adjusts the radii to create a new model.
                        δ = τ3 * δold
                        Δ = max( τ3 * Δold, Δinit )
                        nΓ += 1

                    end

                    # Constructs a new model and set the new relative bounds 'ao' and 'bo'.
                    construct_new_model!( func_list, imin_idx, δ, fi_x, x, a, b, ao, bo, model )                   

                else

                    # Updates the model with TRSBOX information
                    update_model!( t, fi_x, x, model, trsbox_step = true)

                end

            else

                # Updates the model with ALTMOV information
                update_model!( t, fi_x, x, model)

            end

            # Increases iteration counter
            nit += 1

        end

        #---------------------- Preparations to finish execution  ----------------------

        if ( it_flag == :nonspecified ) || ( it_flag == :criticality ) || ( exit_flag == :numerical_error )

            verif_new_info = false

        end

        # Verifies the information obtained in the last iteration.
        if verif_new_info && ( fi_x < model.fval[ model.kopt[] + 1 ] )

            # Calculates the value of fmin if there is enough computational budget and it has not already been calculated. 
            if ( ( maxfun - nf ) ≥ ( r - 1 ) ) && full_calc == false

                fi_x, imin_idx = fmin_partial_eval!( func_list, r, model.imin[], fi_x, x, imin_set )
                nf += r - 1
                full_calc = true

                # Verifies if 'model.imin' belongs to the 'imin_set' of the new point.
                if !( imin_set[ model.imin[] ] )

                    model.imin[] = imin_idx

                end

                # Updates 'model.xopt' and its function value.
                @. model.xopt = x
                model.fval[ model.kopt[] + 1 ] = fi_x

            else

                # Updates 'model.xopt' and its function value.
                @. model.xopt = x
                model.fval[ model.kopt[] + 1 ] = fi_x

            end
            
        end
        
        # Creates the LOWDEROutput.
        output = create_output(model, exit_flag, full_calc, nit, nf, δold, Δold, π)

        # Prints information about the iteration, exit flag and LOWDEROutput.
        print_info(model, output, it_flag, status_flag, full_calc, verbose, δold, Δold, π, ρ, pred_red, real_red, d)

        if saveinfo
    
            save_info!(model, it_flag, status_flag, full_calc, nit, nf, δold, Δold, π, ρ, pred_red, real_red, d, data_file)
            println(data_file, exit_flag)
            close(data_file)
        
        end
        
        return output

    end

end
