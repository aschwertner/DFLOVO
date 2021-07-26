#-------------------------------------------------------------------------------
# Test suite for 'lovo_utils.jl' file
#-------------------------------------------------------------------------------

@testset "lovo_utils" begin
    
    @testset "fmin_eval" begin

        f(x) = ( x[1] + x[2] ) ^ 2.0
        g(x) = x[1] - x[2]
        h(x) = 1.0 - x[1] - x[2]
        l(x) = sin( x[1] + x[2] )
        func_list = [f, g, h, l]
        r = 4
        
        # Test parameters
        y = [0.0, 0.0]
        # Expected return
        fmin_return = 0.0
        imin_return = 1
        # Test
        fmin, imin = fmin_eval(func_list, r, y)
        @test( fmin == fmin_return )
        @test( imin == imin_return )

        # Test parameters
        y = [0.0, 1.0]
        # Expected return
        fmin_return = -1.0
        imin_return = 2
        # Test
        fmin, imin = fmin_eval(func_list, r, y)
        @test( fmin == fmin_return )
        @test( imin == imin_return )

        # Test parameters
        y = [1.0, 1.0]
        # Expected return
        fmin_return = -1.0
        imin_return = 3
        # Test
        fmin, imin = fmin_eval(func_list, r, y)
        @test( fmin == fmin_return )
        @test( imin == imin_return )

        # Test parameters
        y = [0.0, - 1.0]
        # Expected return
        fmin_return = sin(- 1.0)
        imin_return = 4
        # Test
        fmin, imin = fmin_eval(func_list, r, y)
        @test( fmin == fmin_return )
        @test( imin == imin_return )

    end

    @testset "fi_eval" begin

        f(x) = ( x[1] + x[2] ) ^ 2.0
        g(x) = x[1] - x[2]
        h(x) = 1.0 - x[1] - x[2]
        l(x) = sin( x[1] + x[2] )
        func_list = [f, g, h, l]

        # Test parameters
        index = 1
        y = [0.0, - 1.0]
        # Expected return
        value_return = 1.0
        # Test
        value = fi_eval(func_list, index, y)
        @test(value == value_return)
        
        # Test parameters
        index = 2
        y = [1.0, 1.0]
        # Expected return
        value_return = 0.0
        # Test
        value = fi_eval(func_list, index, y)
        @test(value == value_return)

        # Test parameters
        index = 3
        y = [0.0, 1.0]
        # Expected return
        value_return = 0.0
        # Test
        value = fi_eval(func_list, index, y)
        @test(value == value_return)

        # Test parameters
        index = 4
        y = [0.0, 0.0]
        # Expected return
        value_return = 0.0
        # Test
        value = fi_eval(func_list, index, y)
        @test(value == value_return)

    end 

    @testset "verify_j_imin" begin

        f(x) = ( x[1] + x[2] ) ^ 2.0
        g(x) = x[1] - x[2]
        h(x) = 1.0 - x[1] - x[2]
        l(x) = sin( x[1] + x[2] )
        func_list = [f, g, h, l]

        # Test parameters
        y = [0.0, 0.0]
        fmin_y = 0.0
        # Expected return
        test_1_return = true
        test_2_return = true
        test_3_return = false
        test_4_return = true
        # Test
        @test( verify_j_imin(func_list, 1, y, fmin_y) == test_1_return )
        @test( verify_j_imin(func_list, 2, y, fmin_y) == test_2_return )
        @test( verify_j_imin(func_list, 3, y, fmin_y) == test_3_return )
        @test( verify_j_imin(func_list, 4, y, fmin_y) == test_4_return )

        # Test parameters
        y = [0.0, 1.0]
        fmin_y = - 1.0
        # Expected return
        test_1_return = false
        test_2_return = true
        test_3_return = false
        test_4_return = false
        # Test
        @test( verify_j_imin(func_list, 1, y, fmin_y) == test_1_return )
        @test( verify_j_imin(func_list, 2, y, fmin_y) == test_2_return )
        @test( verify_j_imin(func_list, 3, y, fmin_y) == test_3_return )
        @test( verify_j_imin(func_list, 4, y, fmin_y) == test_4_return )

        # Test parameters
        y = [1.0, 1.0]
        fmin_y = - 1.0
        # Expected return
        test_1_return = false
        test_2_return = false
        test_3_return = true
        test_4_return = false
        # Test
        @test( verify_j_imin(func_list, 1, y, fmin_y) == test_1_return )
        @test( verify_j_imin(func_list, 2, y, fmin_y) == test_2_return )
        @test( verify_j_imin(func_list, 3, y, fmin_y) == test_3_return )
        @test( verify_j_imin(func_list, 4, y, fmin_y) == test_4_return )

        # Test parameters
        y = [0.0, - 1.0]
        fmin_y = sin(- 1.0)
        # Expected return
        test_1_return = false
        test_2_return = false
        test_3_return = false
        test_4_return = true
        # Test
        @test( verify_j_imin(func_list, 1, y, fmin_y) == test_1_return )
        @test( verify_j_imin(func_list, 2, y, fmin_y) == test_2_return )
        @test( verify_j_imin(func_list, 3, y, fmin_y) == test_3_return )
        @test( verify_j_imin(func_list, 4, y, fmin_y) == test_4_return )

    end

end