#-------------------------------------------------------------------------------
# Test suite for 'lowder' function
#-------------------------------------------------------------------------------

@testset "lowder" begin

    @testset "linear models" begin
        
        # Linear 01
        function f(x)
            return 2.0 * x[1] - 5.0 * x[2]
        end

        # Linear 02
        function g(x)
            return - 1.5 * x[1] + x[2]
        end

        fmin_list = [f, g]
        x = [3.0, 3.0]
        a = [0.0, 0.0]
        b = [5.0, 5.0]
        Δ = 1.5
        δ = 1.0

        sol = LOWDER.lowder(fmin_list, x, a, b, δ, Δ, m = 3)

        @test(sol.status == :success)
        @test(sol.solution == [0.0, 5.0])
        
    end

end

