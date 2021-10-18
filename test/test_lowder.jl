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

        sol = LOWDER.lowder(fmin_list, x, a, b, δ, Δ; m = 3)

        @test( sol.status == :success )
        @test( sol.index == 1 )
        @test( sol.solution == [0.0, 5.0] )

        # --------------------------------------------------------------------------------
        # Test 02 - Two linear functions (sol.index == 2)
        # --------------------------------------------------------------------------------

        x = [4.0, 1.0]

        sol = LOWDER.lowder(fmin_list, x, a, b, δ, Δ; m = 3)

        @test( sol.status == :success )
        @test( sol.index == 2 )
        @test( sol.solution == [5.0, 0.0] )

        # --------------------------------------------------------------------------------
        # Test 03 - A linear function and a quadratic function (sol.index == 1)
        # --------------------------------------------------------------------------------

        function f3(x)
            return - 0.5 * x[1] + 10 * x[2]
        end
        
        function f4(x)
            return x[1] ^ 2.0 + x[2] ^ 2.0
        end

        fmin_list = [f3, f4]
        a = [0.0, 0.0]
        b = [5.0, 5.0]
        Δ = 1.2
        δ = 1.0

        x = [4.0, 1.0]

        sol = LOWDER.lowder(fmin_list, x, a, b, δ, Δ; m = 3)

        @test( sol.status == :success )
        @test( sol.index == 1 )
        @test( sol.solution == [5.0, 0.0] )

        # --------------------------------------------------------------------------------
        # Test 04 - A linear function and a quadratic function (sol.index == 2)
        # --------------------------------------------------------------------------------

        x = [1.0, 1.0]

        sol = LOWDER.lowder(fmin_list, x, a, b, δ, Δ; m = 3)

        @test( sol.status == :success )
        @test( sol.index == 2 )
        @test( sol.solution == [0.0, 0.0] )

        # --------------------------------------------------------------------------------
        # Test 05 - Two quadratic functions (sol.index == 1)
        # --------------------------------------------------------------------------------

        function f5(x)
            return - x[2] ^ 2.0 + x[1]
        end
        
        function f6(x)
            return 0.1 * x[1] ^ 2.0 - x[1] * x[2]
        end

        fmin_list = [f5, f6]
        a = [0.0, 0.0]
        b = [5.0, 5.0]
        Δ = 1.0
        δ = 0.5

        x = [1.0, 3.0]

        sol = LOWDER.lowder(fmin_list, x, a, b, δ, Δ; m = 3)

        @test( sol.status == :success )
        @test( sol.index == 1 )
        @test( sol.solution == [0.0, 5.0] )

        # --------------------------------------------------------------------------------
        # Test 06 - Two quadratic functions (sol.index == 2)
        # --------------------------------------------------------------------------------

        x = [5.0, 1.0]

        sol = LOWDER.lowder(fmin_list, x, a, b, δ, Δ; m = 3)

        @test( sol.status == :success )
        @test( sol.index == 2 )
        @test( sol.solution == [5.0, 5.0] )

        # --------------------------------------------------------------------------------
        # Test 07 - Rosenbrock & Freudenstein and Roth functions 
        #           (full_active_set and η = 0.0)
        # --------------------------------------------------------------------------------

        function f7(x)
            return ( 10.0 * ( x[2] - x[1] ^ 2.0 ) ) ^ 2.0 + ( 1.0 - x[1] ) ^ 2.0
        end
        
        function f8(x)
            return ( -13.0 + x[1] + ( ( 5.0 - x[2] ) * x[2] - 2.0 ) * x[2] ) ^ 2.0 + ( -29.0 + x[1] + ( ( x[2] + 1.0 ) * x[2] - 14.0 ) * x[2] ) ^ 2.0
        end

        fmin_list = [f7, f8]
        a = [0.0, 0.0]
        b = [5.0, 5.0]
        Δ = 1.2
        δ = 1.0

        x = [5.0, 0.0]

        sol = LOWDER.lowder(fmin_list, x, a, b, δ, Δ; m = 3, η = 0.0)

        @test( sol.status == :success )
        @test( sol.index == 2 )
        @test( sol.solution == [5.0, 0.0] )

        # --------------------------------------------------------------------------------
        # Test 08 - Rosenbrock & Freudenstein and Roth functions (full_active_set)
        # --------------------------------------------------------------------------------

        x = [5.0, 0.0]

        sol = LOWDER.lowder(fmin_list, x, a, b, δ, Δ; m = 3)

        @test( sol.status == :success )
        @test( sol.index == 2 )
        @test( sol.solution == [5.0, 0.0] )

        # --------------------------------------------------------------------------------
        # Test 09 - Rosenbrock & Freudenstein and Roth functions
        # --------------------------------------------------------------------------------

        x = [5.0, 3.0]
        Δ = 2.0
        δ = 1.0

        sol = LOWDER.lowder(fmin_list, x, a, b, δ, Δ; m = 3)

        @test( sol.status == :success )
        @test( sol.index == 2 )
        @test( sol.solution == [5.0, 4.0] )

        # --------------------------------------------------------------------------------
        # Test 10 - Rosenbrock & Freudenstein and Roth functions
        # --------------------------------------------------------------------------------

        x = [5.0, 5.0]
        Δ = 2.0
        δ = 1.0

        sol = LOWDER.lowder(fmin_list, x, a, b, δ, Δ; m = 3)

        @test( sol.status == :success )
        @test( sol.index == 2 )
        @test( sol.solution == [5.0, 4.0] )

        # --------------------------------------------------------------------------------
        # Test 11 - Linear and Quadratic functions
        # --------------------------------------------------------------------------------

        function f9(x)
            return x[1] ^ 2.0 + x[2] ^ 2.0
        end

        function f10(x)
            return 0.5 * x[1] + 0.1 * x[2] + 1.0
        end

        fmin_list = [f9, f10]
        a = [- 1.0, - 2.0]
        b = [1.5, 2.0]
        Δ = 0.5
        δ = 0.5

        x = [1.5, 1.0]

        sol = LOWDER.lowder(fmin_list, x, a, b, δ, Δ; m = 3)

        @test( sol.status == :success )
        @test( sol.index == 1 )
        @test( sol.solution == [0.0, 0.0] )

    end

end

