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

        # Trust-region
        Δ = 10.0

        # --------------------------------------------------------------------------------
        # Test 01 - Large trust-region and full active set
        # --------------------------------------------------------------------------------
        
        # Test specifications
        model.xopt .= [0.0, 0.0]
        model.fval[1] = model.c[] + dot( model.g, model.xopt )  # Interpolation model

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

        active_set = zeros(Bool, 2)
        x = zeros(Float64, 2)
        d = zeros(Float64, 2)
        v_aux = zeros(Float64, 2)

        status = LOWDER.trsbox!(model, Δ, l, u, active_set, x, d, v_aux)
        @test( x == [0.0, 0.0] )
        @test( d == [0.0, - 4.0] )
        @test( active_set == [true, true] )
        @test( status == :full_active_set )

        # --------------------------------------------------------------------------------
        # Test 04 - Large trust-region and empty active set
        # --------------------------------------------------------------------------------

        # Test specifications
        model.xopt .= [2.0, 2.0]
        model.fval[1] = model.c[] + dot( model.g, model.xopt )

        active_set = zeros(Bool, 2)
        x = zeros(Float64, 2)
        d = zeros(Float64, 2)
        v_aux = zeros(Float64, 2)

        status = LOWDER.trsbox!(model, Δ, l, u, active_set, x, d, v_aux)
        @test( x == [0.0, 0.0] )
        @test( d == [- 2.0, - 2.0] )
        @test( active_set == [true, true] )
        @test( status == :full_active_set )

        # --------------------------------------------------------------------------------
        # Test 05 - Large trust-region and empty active set
        # --------------------------------------------------------------------------------

        # Test specifications
        model.xopt .= [4.0, 4.0]
        model.fval[1] = model.c[] + dot( model.g, model.xopt )

        active_set = zeros(Bool, 2)
        x = zeros(Float64, 2)
        d = zeros(Float64, 2)
        v_aux = zeros(Float64, 2)

        status = LOWDER.trsbox!(model, Δ, l, u, active_set, x, d, v_aux)
        @test( x == [0.0, 0.0] )
        @test( d == [- 4.0, - 4.0] )
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
        l = [- 2.0, - 2.0]
        u = [6.0, 6.0]

        # Trust-region
        Δ = 1.0

        # --------------------------------------------------------------------------------
        # Test 06 - Small trust-region and empty active set (full cauchy step)
        # --------------------------------------------------------------------------------

        # Test specifications
        model.xopt .= [0.0, 0.0]
        model.fval[2] = model.c[] + dot( model.g, model.xopt )

        active_set = zeros(Bool, 2)
        x = zeros(Float64, 2)
        d = zeros(Float64, 2)
        v_aux = zeros(Float64, 2)

        status = LOWDER.trsbox!(model, Δ, l, u, active_set, x, d, v_aux)
        @test( isapprox(x, (- Δ / norm(model.g)) * model.g) )
        @test( isapprox(d, (- Δ / norm(model.g)) * model.g) )
        @test( active_set == [false, false] )
        @test( status == :trust_region_boundary )

        # --------------------------------------------------------------------------------
        # Test 07 - Small trust-region and full active set
        # --------------------------------------------------------------------------------

        # Test specifications
        model.xopt .= [- 2.0, - 2.0]
        model.fval[2] = model.c[] + dot( model.g, model.xopt )

        active_set = zeros(Bool, 2)
        x = zeros(Float64, 2)
        d = zeros(Float64, 2)
        v_aux = zeros(Float64, 2)

        status = LOWDER.trsbox!(model, Δ, l, u, active_set, x, d, v_aux)
        @test( x == [- 2.0, - 2.0] )
        @test( d == [0.0, 0.0] )
        @test( active_set == [true, true] )
        @test( status == :full_active_set )

        # --------------------------------------------------------------------------------
        # Test 08 - Small trust-region and partial active set
        # --------------------------------------------------------------------------------

        # Test specifications
        model.xopt .= [- 2.0, 3.0]
        model.fval[2] = model.c[] + dot( model.g, model.xopt )

        active_set = zeros(Bool, 2)
        x = zeros(Float64, 2)
        d = zeros(Float64, 2)
        v_aux = zeros(Float64, 2)

        status = LOWDER.trsbox!(model, Δ, l, u, active_set, x, d, v_aux)
        @test( x == [- 2.0, 2.0] )
        @test( d == [0.0, - 1.0] )
        @test( active_set == [true, false] )
        @test( status == :trust_region_boundary )

        # --------------------------------------------------------------------------------
        # Test 09 - Small trust-region and partial active set
        # --------------------------------------------------------------------------------

        # Test specifications
        model.xopt .= [2.0, -2.0]
        model.fval[2] = model.c[] + dot( model.g, model.xopt )

        active_set = zeros(Bool, 2)
        x = zeros(Float64, 2)
        d = zeros(Float64, 2)
        v_aux = zeros(Float64, 2)

        status = LOWDER.trsbox!(model, Δ, l, u, active_set, x, d, v_aux)
        @test( x == [1.0, - 2.0] )
        @test( d == [- 1.0, 0.0] )
        @test( active_set == [false, true] )
        @test( status == :trust_region_boundary )

        # --------------------------------------------------------------------------------
        # Test 10 - Small trust-region and partial active set
        # --------------------------------------------------------------------------------

        # Test specifications
        model.xopt .= [6.0, - 1.0]
        model.fval[2] = model.c[] + dot( model.g, model.xopt )

        active_set = zeros(Bool, 2)
        x = zeros(Float64, 2)
        d = zeros(Float64, 2)
        v_aux = zeros(Float64, 2)

        status = LOWDER.trsbox!(model, Δ, l, u, active_set, x, d, v_aux)
        @test( isapprox(x, model.xopt + (- Δ / norm(model.g)) * model.g) )
        @test( isapprox(d, (- Δ / norm(model.g)) * model.g) )
        @test( active_set == [false, false] )
        @test( status == :trust_region_boundary )

        # --------------------------------------------------------------------------------
        # --------------------------------------------------------------------------------
        # --------------------------------------------------------------------------------

        # Model specifications
        model = LinearModel(2)
        model.c[] = 0.0
        model.g .= [1.0, 2.0]
        model.kopt[] = 1

        # Box constraints
        l = [- 5.0, - 3.0]
        u = [6.0, 10.0]

        # Trust-region
        Δ = 4.0

        # --------------------------------------------------------------------------------
        # Test 11 - Intersection between the trust-region and 
        #           the lower limits of the y-axis of the box.
        # --------------------------------------------------------------------------------

        # Test specifications
        model.xopt .= [0.0, 0.0]
        model.fval[2] = model.c[] + dot( model.g, model.xopt )

        active_set = zeros(Bool, 2)
        x = zeros(Float64, 2)
        d = zeros(Float64, 2)
        v_aux = zeros(Float64, 2)

        status = LOWDER.trsbox!(model, Δ, l, u, active_set, x, d, v_aux)
        @test( x == [- sqrt(7.0), - 3.0] )
        @test( d == [- sqrt(7.0), - 3.0] )
        @test( active_set == [false, true] )
        @test( status == :trust_region_boundary )

        # --------------------------------------------------------------------------------
        # Test 12 - Intersection between the trust-region and 
        #           the lower limits of the x and y-axis of the box.
        # --------------------------------------------------------------------------------

        # Test specifications
        model.xopt .= [- 2.0, - 2.0]
        model.fval[2] = model.c[] + dot( model.g, model.xopt )

        active_set = zeros(Bool, 2)
        x = zeros(Float64, 2)
        d = zeros(Float64, 2)
        v_aux = zeros(Float64, 2)

        status = LOWDER.trsbox!(model, Δ, l, u, active_set, x, d, v_aux)
        @test( x == [- 5.0, - 3.0] )
        @test( d == [- 3.0, - 1.0] )
        @test( active_set == [true, true] )
        @test( status == :full_active_set )

        # --------------------------------------------------------------------------------
        # Test 13 - Intersection between the trust-region and 
        #           the lower limits of the x-axis of the box.
        # --------------------------------------------------------------------------------

        # Test specifications
        model.xopt .= [-4.0, 5.0]
        model.fval[2] = model.c[] + dot( model.g, model.xopt )

        active_set = zeros(Bool, 2)
        x = zeros(Float64, 2)
        d = zeros(Float64, 2)
        v_aux = zeros(Float64, 2)

        status = LOWDER.trsbox!(model, Δ, l, u, active_set, x, d, v_aux)
        @test( x == [- 5.0, 5.0 - sqrt(15)] )
        @test( d == [- 1.0, - sqrt(15)] )
        @test( active_set == [true, false] )
        @test( status == :trust_region_boundary )

        # --------------------------------------------------------------------------------
        # Test 14 - Intersection between the trust-region and 
        #           the lower limits of the x and y-axis of the box (short step).
        # --------------------------------------------------------------------------------

        # Test specifications
        model.xopt .= [- 4.0, - 2.0]
        model.fval[2] = model.c[] + dot( model.g, model.xopt )

        active_set = zeros(Bool, 2)
        x = zeros(Float64, 2)
        d = zeros(Float64, 2)
        v_aux = zeros(Float64, 2)

        status = LOWDER.trsbox!(model, Δ, l, u, active_set, x, d, v_aux)
        @test( x == [- 5.0, - 3.0] )
        @test( d == [- 1.0, - 1.0] )
        @test( active_set == [true, true] )
        @test( status == :full_active_set )

        # --------------------------------------------------------------------------------
        # Test 15 - Full active set
        # --------------------------------------------------------------------------------

        # Test specifications
        model.xopt .= [- 5.0, - 3.0]
        model.fval[2] = model.c[] + dot( model.g, model.xopt )

        active_set = zeros(Bool, 2)
        x = zeros(Float64, 2)
        d = zeros(Float64, 2)
        v_aux = zeros(Float64, 2)

        status = LOWDER.trsbox!(model, Δ, l, u, active_set, x, d, v_aux)
        @test( x == [- 5.0, - 3.0] )
        @test( d == [ 0.0, 0.0] )
        @test( active_set == [true, true] )
        @test( status == :full_active_set )

        # --------------------------------------------------------------------------------
        # --------------------------------------------------------------------------------
        # --------------------------------------------------------------------------------

        # Model specifications
        model = LinearModel(2)
        model.c[] = 970.0
        model.g .= [-66.0, 692.0]
        model.kopt[] = 0

        # Box constraints
        l = [0.0, 0.0]
        u = [5.0, 5.0]

        # Trust-region
        Δ = 1.2

        # --------------------------------------------------------------------------------
        # Test 16 - Full active set
        # --------------------------------------------------------------------------------

        # Test specifications
        model.xopt .= [5.0, 0.0]
        model.fval[model.kopt[] + 1] = model.c[] + dot( model.g, model.xopt )

        active_set = zeros(Bool, 2)
        x = zeros(Float64, 2)
        d = zeros(Float64, 2)
        v_aux = zeros(Float64, 2)

        status = LOWDER.trsbox!(model, Δ, l, u, active_set, x, d, v_aux)
        @test( x == [ 5.0, 0.0 ] )
        @test( d == [ 0.0, 0.0 ] )
        @test( active_set == [true, true] )
        @test( status == :full_active_set )

    end

    @testset "altmov" begin

        # --------------------------------------------------------------------------------
        # Test 01
        # --------------------------------------------------------------------------------

        # Model specifications
        model = LinearModel(2)
        model.c[] = 970.0
        model.g .= [1.0, 2.0]
        model.kopt[] = 0
        model.xbase .= zeros(Float64, 2)
        model.xopt .= zeros(Float64, 2)
        model.Y .= Matrix{Float64}(I, 2, 2)
        model.dst .= [1.0, 1.0]
        model.factorsY .= Matrix{Float64}(I, 2, 2)
        qrY = qr!(model.factorsY, Val(true))
        model.τY .= qrY.τ
        model.jpvtY .= qrY.jpvt

        # Box constraints
        l = [-5.0, 10.0]
        u = [20.0, 15.0]

        # Sample set radius
        δ = 2.0

        # Test specifications
        idx_t = 1

        active_set = zeros(Bool, 2)
        x = zeros(Float64, 2)
        d = zeros(Float64, 2)
        v = zeros(Float64, 2)
        w = zeros(Float64, 2)

        status = LOWDER.altmov!(model, idx_t, δ, l, u, x, d, v, w, active_set)
        @test( x == [2.0, 0.0] )
        @test( d == [2.0, 0.0] )
        @test( status == :usual_altmov )

        # --------------------------------------------------------------------------------
        # Test 02
        # --------------------------------------------------------------------------------

        # Model specifications
        model = LinearModel(2)
        model.c[] = 0.0
        model.g .= [-66.0, 692.0]
        model.kopt[] = 0
        model.xbase .= [5.0, 0.0]
        model.xopt .= [5.0, 0.0]
        model.Y .= [-1.0 0.0; 0.0 1.0]
        model.dst .= [1.0, 1.0]
        model.factorsY .= [-1.0 0.0; 0.0 1.0]
        qrY = qr!(model.factorsY, Val(true))
        model.τY .= qrY.τ
        model.jpvtY .= qrY.jpvt

        # Box constraints
        l = [ 0.0, 0.0 ]
        u = [ 5.0, 5.0]

        # Sample set radius
        δ = 1.0

        # Test specifications
        idx_t = 1

        active_set = zeros(Bool, 2)
        x = zeros(Float64, 2)
        d = zeros(Float64, 2)
        v = zeros(Float64, 2)
        w = zeros(Float64, 2)

        status = LOWDER.altmov!(model, idx_t, δ, l, u, x, d, v, w, active_set)
        @test( x == [4.0, 0.0] )
        @test( d == [-1.0, 0.0] )
        @test( status == :usual_altmov )

    end

end
