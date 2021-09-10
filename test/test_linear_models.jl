@testset "Linear models" begin

    @testset "trsbox" begin

        model = LinearModel(2)

        model.c[] = 0.0
        model.g .= [1.0, 2.0]
        model.kopt[] = 1
        model.xbase .= [0.0, 0.0]
        model.Y[1, :] .= [3.0, 1.0]

        l = [0.0, 0.0]
        u = [4.0, 4.0]

        # Large TR

        Δ = 10.0

        active_set = zeros(Bool, 4)
        x = zeros(Float64, 2)
        d = zeros(Float64, 2)

        trsbox!(model, Δ, l, u, active_set, x, d)

        @test(x == [0.0, 0.0])
        @test(d == [-3.0, -1.0])
        @test(active_set == [true, true, false, false])

        # Small TR

        active_set = zeros(Bool, 4)
        x = zeros(Float64, 2)
        d = zeros(Float64, 2)

        Δ = 1.0
        
        trsbox!(model, Δ, l, u, active_set, x, d)

        dsol = - model.g / norm(model.g) 
        
        @test(x == model.xopt - dsol)
        @test(d == dsol)
        @test(active_set == [false, false, false, false])

    end

end
