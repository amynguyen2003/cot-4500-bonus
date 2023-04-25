import numpy as np

#Givens for 1 and 2 
tol = 1e-6
iters = 50
matrix = np.array([[3, 1, 1], [1, 4, 1], [2, 3, 7]])
column = np.array([1, 3, 0])

#1
def func(x, x_0, tol):
    return (max(abs(x - x_0))) / (max(abs(x_0)) + tol)

def gau_sei(matrix, column, tol, iterations):
    length = len(column)
    x = np.zeros((length), dtype=np.double)
    k = 1

    while (k <= iters):
        x_0 = x.copy()
        for i in range(length):

            sum_one = sum_two = 0

            for j in range(i):
                sum_one += (matrix[i][j] * x[j])

            for j in range(i + 1, length):
                sum_two += (matrix[i][j] * (x_0[j]))

            x[i] = (1 / matrix[i][i]) * (-sum_one - sum_two + column[i])

            if (func(x, x_0, tol) < tol):
                return k

        k = k + 1

    return k

# 2
def jacobi(matrix, column, tol, iters):
    length = len(column)
    x = np.zeros((length), dtype=np.double)
    k = 1

    while (k <= iters):
        x_0 = x.copy()
        for i in range(length):

            sum = 0

            for j in range(length):

                if j != i:
                    sum += (matrix[i][j] * x_0[j])

            x[i] = (1 / matrix[i][i]) * (-sum + column[i])

            if (func(x, x_0, tol) < tol):
                return k

        k = k + 1

    return k

# 3
guess: float = 0.5
toler: float = .000001
seq: str = "x**3 - (x**2) + 2"

def deriv(value):
    return (3 * value * value) - (2 * value)

def newt_raph(guess: float, toler: float, seq: str):
    count = 0
    x = guess
    f = eval(seq)
    derivative = deriv(guess)
    approx: float = f / derivative

    while(abs(approx) >= toler):
        x = guess
        f = eval(seq)
        derivative = deriv(guess)
        approx = f / derivative
        guess -= approx
        count += 1
    return count

# 4
def div_diff_method(matrix: np.array):
    for i in range(2, len(matrix)):
        for j in range(2, i + 2):
            
            if j >= len(matrix[i]) or matrix[i][j] != 0:
                continue
            
            numerator = matrix[i][j - 1] - matrix[i - 1][j - 1]
            denominator = matrix[i][0] - matrix[i - j + 1][0]
            fraction = numerator / denominator
            matrix[i][j] = fraction
    return matrix

def hermite_interpolation():
    x_pts = [0, 1, 2]
    y_pts = [1, 2, 4]
    f_prime = [1.06, 1.23, 1.55]
    size = len(x_pts) 
    matrix = np.zeros((size * 2, size * 2))
    index = 0
    
    for x in range(0, size * 2, 2):
        matrix[x][0] = x_pts[index]
        matrix[x + 1][0] = x_pts[index]
        index = index + 1
        
    index = 0
    for y in range(0, size * 2, 2):
        matrix[y][1] = y_pts[index]
        matrix[y + 1][1] = y_pts[index]
        index = index + 1

    index = 0
    for i in range(1, size * 2, 2):
        matrix[i][2] = f_prime[index]
        index = index + 1
    
    filled_matrix = div_diff_method(matrix)

    print(filled_matrix)

hermite_interpolation()

# 5
in_pt = 0.5
a = 0
b = 3
n = 100  

def function(t, y):
    return y - (t**3)

def modified_eulers(in_pt, a, b, n):
    h = (b - a) / n
    t = a
    w = in_pt
    
    for i in range(n):
        w = w + ((h / 2) * (function(t, w) + 
                            function(t + h, w + 
                            (h * function(t, w)))))
        t = t + h
    
    return w

#Printing results

print(gau_sei(matrix, column, tol, iters))
print("\n")
print(jacobi(matrix, column, tol, iters))
print("\n")
print(newt_raph(guess, toler, seq))
print("\n")
hermite_interpolation()
print("\n")
print("%.5f" %modified_eulers(in_pt, a, b, n))