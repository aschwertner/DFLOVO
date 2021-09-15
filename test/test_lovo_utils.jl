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
        set = zeros(Bool, r)
        # Expected return
        fmin_return = 0.0
        imin_return = 1
        # Test
        fmin, imin = LOWDER.fmin_eval!(func_list, r, y, set)
        @test( fmin == fmin_return )
        @test( imin == imin_return )
        @test( set == [true, true, false, true] )

        # Test parameters
        y = [0.0, 1.0]
        set = zeros(Bool, r)
        # Expected return
        fmin_return = -1.0
        imin_return = 2
        # Test
        fmin, imin = LOWDER.fmin_eval!(func_list, r, y, set)
        @test( fmin == fmin_return )
        @test( imin == imin_return )
        @test( set == [false, true, false, false] )

        # Test parameters
        y = [1.0, 1.0]
        set = zeros(Bool, r)
        # Expected return
        fmin_return = -1.0
        imin_return = 3
        # Test
        fmin, imin = LOWDER.fmin_eval!(func_list, r, y, set)
        @test( fmin == fmin_return )
        @test( imin == imin_return )
        @test( set == [false, false, true, false] )

        # Test parameters
        y = [0.0, - 1.0]
        set = zeros(Bool, r)
        # Expected return
        fmin_return = sin(- 1.0)
        imin_return = 4
        # Test
        fmin, imin = LOWDER.fmin_eval!(func_list, r, y, set)
        @test( fmin == fmin_return )
        @test( imin == imin_return )
        @test( set == [false, false, false, true] )

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
        value = LOWDER.fi_eval(func_list, index, y)
        @test(value == value_return)
        
        # Test parameters
        index = 2
        y = [1.0, 1.0]
        # Expected return
        value_return = 0.0
        # Test
        value = LOWDER.fi_eval(func_list, index, y)
        @test(value == value_return)

        # Test parameters
        index = 3
        y = [0.0, 1.0]
        # Expected return
        value_return = 0.0
        # Test
        value = LOWDER.fi_eval(func_list, index, y)
        @test(value == value_return)

        # Test parameters
        index = 4
        y = [0.0, 0.0]
        # Expected return
        value_return = 0.0
        # Test
        value = LOWDER.fi_eval(func_list, index, y)
        @test(value == value_return)

    end 

    @testset "fmin_partial_eval" begin

        f(x) = ( x[1] + x[2] ) ^ 2.0
        g(x) = x[1] - x[2]
        h(x) = 1.0 - x[1] - x[2]
        l(x) = sin( x[1] + x[2] )
        func_list = [f, g, h, l]
        r = 4
        set = zeros(Bool, r)

        # Test parameters
        y = [0.0, 0.0]
        
        idx = 1
        fi_y = 0.0

        # Test 01
        fmin, imin = LOWDER.fmin_partial_eval!( func_list, r, idx, fi_y, y, set )
        @test( fmin == 0.0 )
        @test( imin == 1 )
        @test( set == [true, true, false, true] )

        # Test parameters
        idx = 2
        fi_y = 0.0

        # Test 02
        fmin, imin = LOWDER.fmin_partial_eval!( func_list, r, idx, fi_y, y, set )
        @test( fmin == 0.0 )
        @test( imin == 2 )
        @test( set == [true, true, false, true] )

        # Test parameters
        idx = 3
        fi_y = 1.0

        # Test 03
        fmin, imin = LOWDER.fmin_partial_eval!( func_list, r, idx, fi_y, y, set )
        @test( fmin == 0.0 )
        @test( imin == 1 )
        @test( set == [true, true, false, true] )

        # Test parameters
        idx = 4
        fi_y = 0.0

        # Test 04
        fmin, imin = LOWDER.fmin_partial_eval!( func_list, r, idx, fi_y, y, set )
        @test( fmin == 0.0 )
        @test( imin == 4 )
        @test( set == [true, true, false, true] )

        # ---        

        # Test parameters
        y = [0.0, 1.0]
        idx = 3
        fi_y = 0.0

        # Test
        fmin, imin = LOWDER.fmin_partial_eval!( func_list, r, idx, fi_y, y, set )
        @test( fmin == -1.0 )
        @test( imin == 2 )
        @test( set == [false, true, false, false] )

        # --- 

        # Test parameters
        y = [1.0, 1.0]
        idx = 1
        fi_y = 4.0

        # Test
        fmin, imin = LOWDER.fmin_partial_eval!( func_list, r, idx, fi_y, y, set )
        @test( fmin == -1.0 )
        @test( imin == 3 )
        @test( set == [false, false, true, false] )

        # ---        

        # Test parameters
        y = [0.0, -1.0]
        idx = 1
        fi_y = 0.0

        # Test
        fmin, imin = LOWDER.fmin_partial_eval!( func_list, r, idx, fi_y, y, set )
        @test( fmin == l(y) )
        @test( imin == 4 )
        @test( set == [false, false, false, true] )

    end

end