#-------------------------------------------------------------------------------
# Set of useful functions related to LOWDEROutput struct.
#-------------------------------------------------------------------------------

function create_output(
                        model::AbstractModel,
                        nit::Int64,
                        nf::Int64,
                        exit_flag::Int64,
                        full_calc::Bool
                        )

    if ( model.kopt[] == 1 ) || ( full_calc )
        return LOWDEROutput(nit, model.imin[], nf, exit_flag, true, model.fval[1], model.xopt)
    else
        return LOWDEROutput(nit, model.imin[], nf, exit_flag, false, model.fval[model.kopt[]], model.xopt)
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