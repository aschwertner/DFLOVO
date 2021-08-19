struct LOWDEROutput
   
    iter       :: Int64
    index      :: Int64
    nf         :: Int64
    status     :: Int64
    true_val   :: Bool
    f          :: Float64
    solution   :: Vector{Float64}

end

###!!! Escrever uma função para criar a estrutura do LOWDEROutput

function show_output(
                        output::LOWDEROutput
                        )
    
    println("** LOWDER Output **")
    println("Status (.status): ", $(output.status))
    println("Solution (.solution): ", $(output.solution))
    println("Number of iterations (.iter): ", $(output.iter))
    println("Number of function evaluations (.nf): ", $(output.nf))
    println("Objective function value (.f): ", $(output.f))
    println("Trusted objective funtion value (.true_val): ", $(output.true_val))
    println("Index of the function (.index): ", $(output.index))

end