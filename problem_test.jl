using LOWDER

# Rosenbrock function
function f(x)

    return ( 10.0 * ( x[2] - x[1] ^ 2.0 ) ) ^ 2.0 + ( 1.0 - x[1] ) ^ 2.0

end

# Freudenstein and Roth function
function g(x)

    return ( -13.0 + x[1] + ( ( 5.0 - x[2] ) * x[2] - 2.0 ) * x[2] ) ^ 2.0 + ( -29.0 + x[1] + ( ( x[2] + 1.0 ) * x[2] - 14.0 ) * x[2] ) ^ 2.0

end

fmin_list = [f, g]
x = [3.0, 3.0]
a = [0.0, 0.0]
b = [5.0, 5.0]

sol = lowder(fmin_list, x, a, b, 1.0, 1.5, m = 3, maxit = 100, maxfun = 75, verbose = 3)