#-------------------------------------------------------------------------------
# Test suite for 'lowder' function
#-------------------------------------------------------------------------------

@testset "lowder" begin

    @testset "linear models" begin
        
        # --------------------------------------------------------------------------------
        # Test 01 - Two linear functions (sol.index == 1)
        # --------------------------------------------------------------------------------

        function f1(x)
            return 2.0 * x[1] - 5.0 * x[2]
        end

        function f2(x)
            return - 1.5 * x[1] + x[2]
        end

        fmin_list = [f1, f2]
        a = [0.0, 0.0]
        b = [5.0, 5.0]
        Δ = 1.5
        δ = 1.0

        x = [1.0, 2.0]

        sol = LOWDER.lowder(fmin_list, x, a, b, δ, Δ, m = 3)

        @test( sol.status == :success )
        @test( sol.index == 1 )
        @test( sol.solution == [0.0, 5.0] )

        # --------------------------------------------------------------------------------
        # Test 02 - Two linear functions (sol.index == 2)
        # --------------------------------------------------------------------------------

        x = [4.0, 1.0]

        sol = LOWDER.lowder(fmin_list, x, a, b, δ, Δ, m = 3)

        @test( sol.status == :success )
        @test( sol.index == 2 )
        @test( sol.solution == [5.0, 0.0] )

        # --------------------------------------------------------------------------------
        # Test 03 - A linear function and a quadratic function (sol.index == 1)
        # --------------------------------------------------------------------------------

        function g1(x)
            return - 0.5 * x[1] + 10 * x[2]
        end
        
        function g2(x)
            return x[1] ^ 2.0 + x[2] ^ 2.0
        end

        fmin_list = [g1, g2]
        a = [0.0, 0.0]
        b = [5.0, 5.0]
        Δ = 1.2
        δ = 1.0

        x = [4.0, 1.0]

        sol = LOWDER.lowder(fmin_list, x, a, b, δ, Δ, m = 3)

        @test( sol.status == :success )
        @test( sol.index == 1 )
        @test( sol.solution == [5.0, 0.0] )

        # --------------------------------------------------------------------------------
        # Test 04 - A linear function and a quadratic function (sol.index == 2)
        # --------------------------------------------------------------------------------

        x = [1.0, 1.0]

        sol = LOWDER.lowder(fmin_list, x, a, b, δ, Δ, m = 3)

        @test( sol.status == :success )
        @test( sol.index == 2 )
        @test( sol.solution == [0.0, 0.0] )

        # --------------------------------------------------------------------------------
        # Test 05 - Two quadratic functions (sol.index == 1)
        # --------------------------------------------------------------------------------

        function h1(x)
            return - x[2] ^ 2.0 + x[1]
        end
        
        function h2(x)
            return 0.1 * x[1] ^ 2.0 - x[1] * x[2]
        end

        fmin_list = [h1, h2]
        a = [0.0, 0.0]
        b = [5.0, 5.0]
        Δ = 1.0
        δ = 0.5

        x = [1.0, 3.0]

        sol = LOWDER.lowder(fmin_list, x, a, b, δ, Δ, m = 3)

        @test( sol.status == :success )
        @test( sol.index == 1 )
        @test( sol.solution == [0.0, 5.0] )

        # --------------------------------------------------------------------------------
        # Test 06 - Two quadratic functions (sol.index == 2)
        # --------------------------------------------------------------------------------

        x = [5.0, 1.0]

        sol = LOWDER.lowder(fmin_list, x, a, b, δ, Δ, m = 3)

        @test( sol.status == :success )
        @test( sol.index == 2 )
        @test( sol.solution == [5.0, 5.0] )

    end

end

