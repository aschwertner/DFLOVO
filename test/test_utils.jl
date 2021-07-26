#-------------------------------------------------------------------------------
# Test suite for 'utils.jl' file
#-------------------------------------------------------------------------------

@testset "utils" begin
    
    @testset "verify_initial_room" begin
        
        # Test parameters
        n = 2
        Δ = 2.0
        a = [1.0, 3.0]
        b = [5.0, 10.0]
        # Test
        @test( verify_initial_room(n, Δ, a, b) == true )

        # Test parameters
        n = 3
        Δ = 1.0
        a = [0.0, 0.0, 0.0]
        b = [1.0, 5.0, 12.0]
        # Test
        @test( verify_initial_room(n, Δ, a, b) == false )

        # Test parameters
        n = 5
        Δ = 3.0
        a = [- 2.0, - 4.0, 0.0, 1.0, 1.5]
        b = [1.0, - 2.0, 1.0, 2.5, 6.0]
        # Test
        @test( verify_initial_room(n, Δ, a, b) == false )

    end

    @testset "correct_guess_bounds!" begin

        # Test parameters
        n = 2
        Δ = 2.0
        a = [1.0, 3.0]
        b = [5.0, 10.0]
        x = [3.0, 5.0]
        ao = [0.0, 0.0]
        bo = [0.0, 0.0]
        # Expected return
        x_return = [3.0, 5.0]
        ao_return = [- 2.0, - 2.0]
        bo_return = [2.0, 5.0]
        # Test
        correct_guess_bounds!(n, Δ, a, b, x, ao, bo)
        @test( x == x_return )
        @test( ao == ao_return)
        @test( bo == bo_return)

        # Test parameters
        n = 2
        Δ = 2.0
        a = [1.0, 3.0]
        b = [5.0, 10.0]
        x = [0.0, 0.0]
        ao = [0.0, 0.0]
        bo = [0.0, 0.0]
        # Expected return
        x_return = [1.0, 3.0]
        ao_return = [0.0, 0.0]
        bo_return = [4.0, 7.0]
        # Test
        correct_guess_bounds!(n, Δ, a, b, x, ao, bo)
        @test( x == x_return )
        @test( ao == ao_return)
        @test( bo == bo_return)

        # Test parameters
        n = 2
        Δ = 2.0
        a = [1.0, 3.0]
        b = [5.0, 10.0]
        x = [10.0, 10.0]
        ao = [0.0, 0.0]
        bo = [0.0, 0.0]
        # Expected return
        x_return = [5.0, 10.0]
        ao_return = [- 4.0, - 7.0]
        bo_return = [0.0, 0.0]
        # Test
        correct_guess_bounds!(n, Δ, a, b, x, ao, bo)
        @test( x == x_return )
        @test( ao == ao_return)
        @test( bo == bo_return)

        # Test parameters
        n = 2
        Δ = 2.0
        a = [1.0, 3.0]
        b = [5.0, 10.0]
        x = [0.0, 20.0]
        ao = [0.0, 0.0]
        bo = [0.0, 0.0]
        # Expected return
        x_return = [1.0, 10.0]
        ao_return = [0.0, - 7.0]
        bo_return = [4.0, 0.0]
        # Test
        correct_guess_bounds!(n, Δ, a, b, x, ao, bo)
        @test( x == x_return )
        @test( ao == ao_return)
        @test( bo == bo_return)

        # Test parameters
        n = 2
        Δ = 2.0
        a = [1.0, 3.0]
        b = [5.0, 10.0]
        x = [2.0, 5.0]
        ao = [0.0, 0.0]
        bo = [0.0, 0.0]
        # Expected return
        x_return = [3.0, 5.0]
        ao_return = [- 2.0, - 2.0]
        bo_return = [2.0, 5.0]
        # Test
        correct_guess_bounds!(n, Δ, a, b, x, ao, bo)
        @test( x == x_return )
        @test( ao == ao_return)
        @test( bo == bo_return)

        # Test parameters
        n = 2
        Δ = 2.0
        a = [1.0, 3.0]
        b = [5.0, 10.0]
        x = [3.0, 4.0]
        ao = [0.0, 0.0]
        bo = [0.0, 0.0]
        # Expected return
        x_return = [3.0, 5.0]
        ao_return = [- 2.0, - 2.0]
        bo_return = [2.0, 5.0]
        # Test
        correct_guess_bounds!(n, Δ, a, b, x, ao, bo)
        @test( x == x_return )
        @test( ao == ao_return)
        @test( bo == bo_return)

        # Test parameters
        n = 2
        Δ = 2.0
        a = [1.0, 3.0]
        b = [5.0, 10.0]
        x = [2.0, 4.0]
        ao = [0.0, 0.0]
        bo = [0.0, 0.0]
        # Expected return
        x_return = [3.0, 5.0]
        ao_return = [- 2.0, - 2.0]
        bo_return = [2.0, 5.0]
        # Test
        correct_guess_bounds!(n, Δ, a, b, x, ao, bo)
        @test( x == x_return )
        @test( ao == ao_return)
        @test( bo == bo_return)

        # Test parameters
        n = 2
        Δ = 2.0
        a = [1.0, 3.0]
        b = [5.0, 10.0]
        x = [4.0, 5.0]
        ao = [0.0, 0.0]
        bo = [0.0, 0.0]
        # Expected return
        x_return = [3.0, 5.0]
        ao_return = [- 2.0, - 2.0]
        bo_return = [2.0, 5.0]
        # Test
        correct_guess_bounds!(n, Δ, a, b, x, ao, bo)
        @test( x == x_return )
        @test( ao == ao_return)
        @test( bo == bo_return)

        # Test parameters
        n = 2
        Δ = 2.0
        a = [1.0, 3.0]
        b = [5.0, 10.0]
        x = [3.0, 9.0]
        ao = [0.0, 0.0]
        bo = [0.0, 0.0]
        # Expected return
        x_return = [3.0, 8.0]
        ao_return = [- 2.0, - 5.0]
        bo_return = [2.0, 2.0]
        # Test
        correct_guess_bounds!(n, Δ, a, b, x, ao, bo)
        @test( x == x_return )
        @test( ao == ao_return)
        @test( bo == bo_return)

        # Test parameters
        n = 2
        Δ = 2.0
        a = [1.0, 3.0]
        b = [5.0, 10.0]
        x = [4.0, 9.0]
        ao = [0.0, 0.0]
        bo = [0.0, 0.0]
        # Expected return
        x_return = [3.0, 8.0]
        ao_return = [- 2.0, - 5.0]
        bo_return = [2.0, 2.0]
        # Test
        correct_guess_bounds!(n, Δ, a, b, x, ao, bo)
        @test( x == x_return )
        @test( ao == ao_return)
        @test( bo == bo_return)

    end

end