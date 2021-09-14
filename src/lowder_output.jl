#-------------------------------------------------------------------------------
# Set of useful functions related to LOWDEROutput struct.
#-------------------------------------------------------------------------------

function create_output(
                        model::AbstractModel,
                        exit_flag::Symbol,
                        full_calc::Bool,
                        nit::Int64,
                        nf::Int64
                        )

    if ( model.kopt[] == 0 ) || ( full_calc )

        return LOWDEROutput(exit_flag, true, nit, nf, model.imin[], model.fval[1], model.xopt)

    else

        return LOWDEROutput(exit_flag, false, nit, nf, model.imin[], model.fval[ model.kopt[] + 1 ], model.xopt)

    end

end

function show_output(
                        output::LOWDEROutput
                        )
    
    println("** LOWDER Output **")
    println("Status (.status): $( output.status )" )
    println("Solution (.solution): $( output.solution )" )
    println("Number of iterations (.iter): $( output.iter )" )
    println("Number of function evaluations (.nf): $( output.nf )" )
    println("Objective function value (.f): $( output.f )" )
    println("Trusted objective funtion value (.true_val): $( output.true_val )" )
    println("Index of the function (.index): $( output.index )" )

end