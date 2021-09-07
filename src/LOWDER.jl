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
    include("lovo_utils.jl")
    include("utils.jl")
    include("linear_models.jl")
    include("lowder_output.jl")

    function lowder(
                    func_list::Array{Function,1},
                    x::Vector{Float64},
                    a::Vector{Float64}, 
                    b::Vector{Float64},
                    δ::Float64,
                    Δ::Float64;
                    m::Int64=(2 * length(x) + 1),
                    maxit::Int64=5000,
                    maxfun::Int64=(1000 * (length(func_list) + m)),
                    nρmax::Int64=3,
                    Γmax::Int64=1,
                    verbose::Int64=0,
                    δmin::Float64=1.0e-8,
                    πmin::Float64=1.0e-8,
                    β::Float64=1.0,
                    τ1::Float64=0.6,
                    τ2::Float64=1.5,
                    τ3::Float64=2.0,
                    η::Float64=0.1,
                    η1::Float64=0.3,
                    η2::Float64=0.6
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
        
        # Sets some useful constants.
        Δinit = Δ

        # Initializes useful variables, vectors, and matrices.
        nρ = 0                      # Auxiliary counter for simplified 'ρ' calculations.
        nΓ = 0                      # Auxiliary counter for radii adjustments phase.
        nit = 0                     # Counts the number of iterations.
        nf = 0                      # Counts the number of 'f_{i}' function evaluations.
        ao = zeros(Float64, n)      # Difference between the lower bounds 'a' and the center of the sample set, given by 'xbase'.
        bo = zeros(Float64, n)      # Difference between the upper bounds 'b' and the center of the sample set, given by 'xbase'.
        d = zeros(Float64, n)       # TRSBOX or ALTMOV direction.
        active_set = zeros(Bool, n) # Set of active constraints.
        imin_set = zeros(Bool, r)   # Set I_{min}(x)
        full_calc = false

        #-------------------- Preparations for the first iteration ---------------------

        # Create an empty structure called 'model', which can hold a linear or quadratic model.
        if n == (m - 1)
            model = create_linear_model(n)
        end

        # Modifies the initial estimate 'x' to be suitable for building the first model. 
        # Modifies 'ao' and 'bo' to store the 'a-x' and 'b-x' differences, respectively.
        correct_guess_bounds!(n, δ, a, b, x, ao, bo)

        # Computes the value of f_{min}('x') and an index 'imin_idx' ∈ I_{min}('x').
        fbase, imin_idx = fmin_eval!(func_list, r, x, imin_set)

        # Updates de function call counter.
        nf += r

        # Creates the initial model, modifying the previous 'model' structure.
        construct_model!(func_list, imin_idx, δ, fbase, x, ao, bo, model)

        # Updates de function call counter.
        nf += n - 1

        # Returns if 'nf' exceeds 'maxfun'.
        if nf ≥ maxfun

            # Creates the LOWDEROutput.
            output = create_output(model, nit, nf, -2, full_calc)

            # Prints information about the iteration, exit flag and LOWDEROutput.
            print_info(model, output, false, verbose, -2, 0, nit, nf, δ, Δ)
            
            return output
            
        end
        
        # Main loop
        while true

            # Saves old information about radii
            δold = δ
            Δold = Δ

            # Computes the stationarity measure π = ||P_{Ω}(xopt - g) - xopt||
            π = stationarity_measure(model, a, b)

            # Verifies if 'δ' and 'π' are less than or equal to 'δmin' and 'πmin', respectively.
            if ( δ ≤ δmin ) && ( π ≤ πmin )

                exit_flag = 1
                break

            end

            if δ > β * π

                #------------------------------ Criticality phase ------------------------------

                # Sets iteration flag
                it_flag = 1

                # Update parameters
                δ *= τ1
                ρ = 0.0

            else

                #------------------------------- Step calculation ------------------------------

                trsbox!(model, Δ, ao, bo, active_set, x, d)

                #------------------------------- Step acceptance -------------------------------

                mnew = model(x)
                diff = model.fval[model.kopt[]] - mnew
                full_calc = false

                if diff < 0

                    it_flag = 0
                    exit_flag = -3
                    break

                else

                    if nρ < nρmax

                        ρ, fi_x, x_idx = relative_reduction(model, func_list, diff, x)
                        nf +=1

                    else

                        full_calc = true
                        ρ, fi_x, x_idx = relative_reduction!(model, func_list, r, diff, x, imin_set)
                        nf += r - 1

                    end

                end

                #------------------------------- Radii updates ---------------------------------

                if ρ < η1

                    δ *= τ1
                    Δ *= τ1

                elseif ( ρ > η2 ) && ( norm(d) == Δ )

                    δ *= τ2
                    Δ *= τ2

                end

            end

            # Chooses the point 't' that must leave the interpolation set. In the case of an ALTMOV call,
            # it also computes the new point 'x' and the direction 'd'.
            if ρ ≥ η

                t = choose_index_trsbox(model, Δ, x)

                it_flag = 2
                nρ = 0

            else

                t = choose_index_altmov(model)

                altmov_flag = altmov!(model, t, ao, bo, x, d, active_set)

                if altmov_flag
                    it_flag = 3
                else
                    it_flag = 4
                end
                nρ += 1

            end

            #------------------------- Verifies output conditions --------------------------
        
            # Verifies if 'nit' exceeds 'maxit', or if 'nf' exceeds 'maxfun'. 
            if nit ≥ maxit

                exit_flag = -1
                break

            elseif nf ≥ maxfun

                exit_flag = -2
                break

            else

                if verbose ≥ 1
                    print_iteration( full_calc, it_flag, nit, nf, model.imin[], δold, Δold, model.fval[model.kopt[]] )
                end

            end

            #--------------------------- Radii adjustments phase ---------------------------

            if ρ ≥ η

                if full_calc

                    if imin_set[ model.imin[] ]

                        x_idx = model.imin[]

                    end

                else

                    fi_x, x_idx = fmin_partial_eval!( func_list, r, x_idx, fi_x, x, imin_set )
                    nf += r - 1

                    if imin_set[ model.imin[] ]

                        x_idx = model.imin[]

                    end

                end

                if x_idx != model.imin[]

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
                    construct_new_model!( func_list, x_idx, δ, fi_x, x, a, b, ao, bo, model )                   

                else

                    # Updates the model with TRSBOX information
                    update_model!( t, fi_x, x, model, trsbox_step = true)

                end

            else

                # Computes the new function value, increases the function evaluation counter, and updates the model with ALTMOV information
                fi_x = fi_eval( func_list, model.imin[], x)
                nf += 1
                update_model!( t, fi_x, x, model)

            end

            # Increases iteration counter
            nit += 1

        end

        #---------------------- Preparations to finish execution  ----------------------

        if fi_x < model.fval[ model.kopt[] ]

            @. model.xopt = x
            model.fval[ model.kopt[] ] = fi_x

        end
        
        # Creates the LOWDEROutput.
        output = create_output(model, nit, nf, exit_flag, full_calc)

        # Prints information about the iteration, exit flag and LOWDEROutput.
        print_info(model, output, false, verbose, exit_flag, it_flag, nit, nf, δ, Δ)
        
        return output

    end

end
