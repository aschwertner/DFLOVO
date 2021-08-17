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
    include("lovo_utils.jl")
    include("utils.jl")
    include("linear_models.jl")
    include("linear_altmov.jl")

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
                    ρmax::Int64=3,
                    Γmax::Int64=1,
                    δmin::Float64=1.0e-8,
                    πmin::Float64=1.0e-8,
                    β::Float64=1.0,
                    τ1::Float64=0.6,
                    τ2::Float64=1.5,
                    τ3::Float64=2.0,
                    η::Float64=0.1,
                    η1::Float64=0.3,
                    η2::Float64=0.6,
                    verbose::Bool=true
                    )

        # Calculates the search space dimension.
        n = length(x)

        # Calculates the number of functions that make up the objective function fmin
        r = length(func_list)

        # Verify algorithm initialization conditions.
        m_min = n + 1
        m_max = convert(Int64, ( n + 1 ) * ( n + 2 ) / 2)
        @assert m_min ≤ m ≤ m_max "The number of interpolation points 'm' must satisfy m ∈ [$(m_min), $(m_max)]."
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
        
        # Sets some useful constants.
        Δinit = Δ
        #nh = convert(Int64, n * ( n + 1 ) / 2)

        # Initializes useful variables, vectors, and matrices.
        countit = 0                 # Counts the number of iterations.
        countf = 0                  # Counts the number of 'f_i' function evaluations.
        countρ = 0                  # Counts the number of simplified 'ρ' calculations for which the descent direction was not accepted.
        Γ = 0                       # Auxiliary counter for Radii adjustments phase.
        xopt = zeros(n)             # Best point so far.
        ao = zeros(n)               # Difference between the lower bounds 'a' and the center of the sample set, given by 'xbase'.
        bo = zeros(n)               # Difference between the upper bounds 'b' and the center of the sample set, given by 'xbase'.
        d_trs = zeros(n)            # TRSBOX descent direction.
        d_alt = zeros(n)            # ALTMOV descent direction.
        fval = zeros(n + 1)         # Set of the function values of the interpolation points.
        Y = zeros(n, n)             # Set of interpolation points, except the origin 'xbase'.

        #-------------------- Preparations for the first iteration ---------------------

        # Modifies the initial estimate 'x' to be suitable for building the first model. 
        # Modifies 'ao' and 'bo' to store the 'a-x' and 'b-x' differences, respectively.
        correct_guess_bounds!(n, δ, a, b, x, ao, bo)

        # Computes the value of f_min at 'x' and an index 'imin' belonging to the set I_min('x').
        fbase, imin = fmin_eval(func_list, r, x)

        # Updates de function call counter.
        countf += r

        # Builds the interpolation set and computes the funtion values.
        fval[1] = fbase
        kopt = construct_initial_set_linear!(func_list, n, imin, δ, fbase, x, bo, xopt, fval, Y)

        # Updates de function call counter.
        countf += n - 1

        # Defines 'kbase' as the position of 'x' in set 'Y', i.e., 'kbase' is set to 1.
        kbase = 1

        # Saves the objective function value at 'x' in
        fsave = fbase

        # Returns if 'countf' exceeds 'maxfun'.
        if countf ≥ maxfun

            it_flag = 0
            exit_flag = -2

            # Prints information about the iteration.
            if verbose
                print_iteration(countit, countf, it_flag, δ, Δ, fsave)
            end

            # Prints information about the exit flag.
            print_info(exit_flag)

            # Prints additional information
            if kopt != kbase
                add_exit_flag = -11
                print_info(add_exit_flag)
            end
            
            return x, fsave, imin, countit, countf, δ, Δ, it_flag, exit_flag, xopt, fval[kopt]
            
        end

        # Constructs the linear Model
        model = LinearModel(n, imin, δ, x, fval, Y)
        
        while true

            # Saves old information
            δold = δ
            Δold = Δ

            # Computes the stationarity measure π = ||P_{Ω}(xopt - g) - xopt||
            π = stationarity_measure(model, xopt, a, b)

            # Verifies if 'δ' and 'π' are less than or equal to 'δmin' and 'πmin', respectively.
            if ( δ ≤ δmin ) && ( π ≤ πmin )
                exit_flag = 1
                break
            end

            if δ > β * π
                #------------------------------ Criticality phase ------------------------------

                it_flag = 1

                # Update parameters
                δ *= τ1
                ρ = 0.0

            else

                #------------------------------- Step acceptance -------------------------------

                dsc_d, s_cal, idx, ρ, fi_new = relative_reduction(model, func_list, r, kopt, countρ, ρmax, xopt, d_trs, x)
                
                if dsc_d

                    it_flag = 2
                    exit_flag = -3

                    # Prints information about the iteration.
                    if verbose
                        print_iteration(countit, countf, it_flag, δ, Δ, fsave)
                    end

                    # Prints information about the exit flag.
                    print_info(exit_flag)

                    # Prints additional information
                    if kopt != kbase
                        add_exit_flag = -11
                        print_info(add_exit_flag)
                    end
            
                    return model.xbase, fsave, model.imin, countit, countf, δ, Δ, it_flag, exit_flag, xopt, fval[kopt]

                end

                if s_cal
                    countf += 1
                else
                    countf += r
                end

                ### VER
                if countf ≥ maxfun

                    it_flag = 0
                    exit_flag = -2
        
                    # Prints information about the iteration.
                    if verbose
                        print_iteration(countit, countf, it_flag, δ, Δ, fsave)
                    end
        
                    # Prints information about the exit flag.
                    print_info(exit_flag)
        
                    # Prints additional information
                    if kopt != kbase
                        add_exit_flag = -11
                        print_info(add_exit_flag)
                    end
                    
                    return model.xbase, fsave, model.imin, countit, countf, δ, Δ, it_flag, exit_flag, xopt, fval[kopt]
                    
                end

                if ρ ≥ η
                    countρ = 0

                    ### Atualizar o ponto e remover o ponto indicado pelo TRSBOX.
                else
                    countρ += 1
                    ### Calcular a nova direção via ALTMOV e remover o ponto indicado.
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

            # Prints information about the iteration.
            if verbose
                print_iteration(countit, countf, it_flag, δold, Δold, fsave)
            end

            #------------------------- Verifies output conditions --------------------------
        
            # Verifies if 'countit' exceeds 'maxit'.
            if countit ≥ maxit
                exit_flag = -1
                break
            end

            # Verifies if 'countf' exceeds 'maxfun'.
            if countf ≥ maxfun
                exit_flag = -2
                break
            end

            #--------------------------- Radii adjustments phase ---------------------------
            
            if ρ ≥ η
                
                if s_cal
                    # Computes fmin(xnew) and an index in Imin(xnew).
                    fmin_new, imin_new = fmin_partial_eval(func_list, r, idx, fi_new, xnew)
                    countf += r - 1
                    if imin_new != model.imin
                        # Verifies if the index 'model.imin' belongs to the set Imin(xnew).
                        i_imin = verify_index_imin(func_list, model.imin, fmin_new, xnew)
                        countf += 1
                        # If it is true, 'imin_new' is set to 'model.imin'.
                        if i_imin
                            imin_new = imin
                        end
                    end
                else
                    # Verifies if the index 'model.imin' belongs to the set Imin(xnew).
                    i_imin = verify_index_imin(func_list, model.imin, fi_new, xnew)

                    # If it is true, 'imin_new' is set to 'model.imin', otherwise, 'imin_new' is set do idx.
                    if i_imin
                        imin_new = model.imin
                    else
                        imin_new = idx
                        fmin_new = fi_new
                    end
                end

                if imin_new != model.imin

                    if ρ ≥ η1

                        Γ = 0

                    end

                    if Γ ≤ Γmax

                        # Adjusts the radii to create a new model.
                        δ = τ3 * δold
                        Δ = max( τ3 * Δold, Δinit )
                        Γ += 1

                    end
                    
                    # Constructs a new linear model
                    model.imin = imin_new
                    model.fval[1] = fmin_new
                    model = model_from_scratch!(model, func_list, δ, b)

                else

                    # Updates the linear model with TRSBOX information
                    update_model!( model, t_trs, fi_new, d_trs)

                end

            else

                # Updates the linear model with ALTMOV information
                update_model!( model, t_alt, fi_new, d_alt)

            end

            # Increases iteration counter
            countit += 1

        end

        #---------------------- Preparations to finish execution  ----------------------

        # Prints information about the exit flag.
        print_info(exit_flag)

        # Prints additional information
        if kopt != kbase
            add_exit_flag = -12
            print_info(add_exit_flag)
        end

        return model.xbase, fsave, model.imin, countit, countf, δ, Δ, it_flag, exit_flag, xopt, model.fval[kopt]

    end

end
