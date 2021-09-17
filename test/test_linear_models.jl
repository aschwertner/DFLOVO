@testset "linear_models" begin

    @testset "trsbox" begin

        # Model specifications
        model = LinearModel(2)
        model.c[] = 0.0
        model.g .= [1.0, 2.0]
        model.kopt[] = 0

        # Box constraints
        l = [0.0, 0.0]
        u = [4.0, 4.0]

        # --------------------------------------------------------------------------------
        # Test 01 - Large trust-region and full active set
        # --------------------------------------------------------------------------------
        
        # Test specifications
        model.xopt .= [0.0, 0.0]
        model.fval[1] = model.c[] + dot( model.g, model.xopt )  # Interpolation model

        # Large TR
        Δ = 10.0

        active_set = zeros(Bool, 2)
        x = zeros(Float64, 2)
        d = zeros(Float64, 2)
        v_aux = zeros(Float64, 2)

        status = LOWDER.trsbox!(model, Δ, l, u, active_set, x, d, v_aux)
        @test( x == [0.0, 0.0] )
        @test( d == [0.0, 0.0] )
        @test( active_set == [true, true] )
        @test( status == :full_active_set )

        # --------------------------------------------------------------------------------
        # Test 02 - Large trust-region and partial active set
        # --------------------------------------------------------------------------------

        # Test specifications
        model.xopt .= [2.0, 0.0]
        model.fval[1] = model.c[] + dot( model.g, model.xopt )

        # Large TR
        Δ = 10.0

        active_set = zeros(Bool, 2)
        x = zeros(Float64, 2)
        d = zeros(Float64, 2)
        v_aux = zeros(Float64, 2)

        status = LOWDER.trsbox!(model, Δ, l, u, active_set, x, d, v_aux)
        @test( x == [0.0, 0.0] )
        @test( d == [-2.0, 0.0] )
        @test( active_set == [true, true] )
        @test( status == :full_active_set )

        # --------------------------------------------------------------------------------
        # Test 03 - Large trust-region and partial active set
        # --------------------------------------------------------------------------------

        # Test specifications
        model.xopt .= [0.0, 4.0]
        model.fval[1] = model.c[] + dot( model.g, model.xopt )

        # Large TR
        Δ = 10.0

        active_set = zeros(Bool, 2)
        x = zeros(Float64, 2)
        d = zeros(Float64, 2)
        v_aux = zeros(Float64, 2)

        status = LOWDER.trsbox!(model, Δ, l, u, active_set, x, d, v_aux)
        @test( x == [0.0, 0.0] )
        @test( d == [0.0, -4.0] )
        @test( active_set == [true, true] )
        @test( status == :full_active_set )

        # --------------------------------------------------------------------------------
        # Test 04 - Large trust-region and empty active set
        # --------------------------------------------------------------------------------

        # Test specifications
        model.xopt .= [2.0, 2.0]
        model.fval[1] = model.c[] + dot( model.g, model.xopt )

        # Large TR
        Δ = 10.0

        active_set = zeros(Bool, 2)
        x = zeros(Float64, 2)
        d = zeros(Float64, 2)
        v_aux = zeros(Float64, 2)

        status = LOWDER.trsbox!(model, Δ, l, u, active_set, x, d, v_aux)
        @test( x == [0.0, 0.0] )
        @test( d == [-2.0, -2.0] )
        @test( active_set == [true, true] )
        @test( status == :full_active_set )

        # --------------------------------------------------------------------------------
        # Test 05 - Large trust-region and empty active set
        # --------------------------------------------------------------------------------

        # Test specifications
        model.xopt .= [4.0, 4.0]
        model.fval[1] = model.c[] + dot( model.g, model.xopt )

        # Large TR
        Δ = 10.0

        active_set = zeros(Bool, 2)
        x = zeros(Float64, 2)
        d = zeros(Float64, 2)
        v_aux = zeros(Float64, 2)

        status = LOWDER.trsbox!(model, Δ, l, u, active_set, x, d, v_aux)
        @test( x == [0.0, 0.0] )
        @test( d == [-4.0, -4.0] )
        @test( active_set == [true, true] )
        @test( status == :full_active_set )

        # --------------------------------------------------------------------------------
        # --------------------------------------------------------------------------------
        # --------------------------------------------------------------------------------

        # Model specifications
        model = LinearModel(2)
        model.c[] = 0.0
        model.g .= [1.0, 2.0]
        model.kopt[] = 1

        # Box constraints
        l = [-2.0, -2.0]
        u = [6.0, 6.0]

        # --------------------------------------------------------------------------------
        # Test 06 - Small trust-region and full active set
        # --------------------------------------------------------------------------------

        # Test specifications
        model.xopt .= [0.0, 0.0]
        model.fval[2] = model.c[] + dot( model.g, model.xopt )  # Interpolation model

        # Small TR
        Δ = 1.0

        active_set = zeros(Bool, 2)
        x = zeros(Float64, 2)
        d = zeros(Float64, 2)
        v_aux = zeros(Float64, 2)

        status = LOWDER.trsbox!(model, Δ, l, u, active_set, x, d, v_aux)
        @test( isapprox(x, (- Δ / norm(model.g)) * model.g) )
        @test( isapprox(d, (- Δ / norm(model.g)) * model.g) )
        @test( active_set == [false, false] )
        @test( status == :trust_region_boundary )


    end

end
