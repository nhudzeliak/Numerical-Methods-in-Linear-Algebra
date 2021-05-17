import math


def scalar_product(vector1, vector2):
    product = 0
    for i in range(len(vector1)):
        product += vector1[i] * vector2[i]
    return product


def vector_norm(vector):
    return math.sqrt(scalar_product(vector, vector))


def divide_vector_by_scalar(vector, scalar):
    result = []
    for coord in vector:
        result.append(coord / scalar)
    return result


def LU_decomposition(a):
    n = len(a)
    u = [[None]*n for i in range(n)]
    l = [[None]*n for i in range(n)]
    for i in range(n):
        for j in range(n):
            if i < j:
                l[i][j] = 0
            elif i == j:
                l[i][j] = 1
            else:
                u[i][j] = 0

    for j in range(n):
        u[0][j] = a[0][j]
    for i in range(1, n):
        l[i][0] = a[i][0] / a[0][0]
    for i in range(1, n):
        for j in range(1, n):
            if i <= j:
                temp_sum = 0
                for k in range(i):
                    temp_sum += l[i][k] * u[k][j]
                u[i][j] = a[i][j] - temp_sum
            if i > j:
                temp_sum = 0
                for k in range(j):
                    temp_sum += l[i][k] * u[k][j]
                l[i][j] = (a[i][j] - temp_sum) / u[j][j]
    return l, u


def solve_Ly_equals_b(l, b):
    y = []
    for i in range(len(b)):
        temp_sum = 0
        for k in range(i):
            temp_sum += l[i][k] * y[k]
        y.append(b[i] - temp_sum)
    return y


def solve_Ux_equals_y(u, y):
    x = [None] * len(y)
    for i in range(len(y)-1, -1, -1):
        temp_sum = 0
        for k in range(i+1, len(y)):
            temp_sum += u[i][k] * x[k]
        try:
            x[i] = (y[i] - temp_sum) / u[i][i]
        except ZeroDivisionError:
            return None
    return x


def inverse_iteration(A, y0, eps, lmbd0=0):
    """Returns iteration approximations for
    minimum eigenvalue and corresponding eigenvector of given matrix
    implementing inverse iteration algorithm
    :param A: input matrix
    :param y0: first approximation of eigenvector
    :param eps: accuracy
    :param lmbd0: first approximation of eigenvalue
    :return: list of tuples with eigenpairs on each iteration
    """
    s0 = scalar_product(y0, y0)
    x0 = divide_vector_by_scalar(y0, vector_norm(y0))
    approx_pairs = [tuple((lmbd0, x0))]
    continue_iter = True
    l, u = LU_decomposition(A)
    while continue_iter:
        z = solve_Ly_equals_b(l, x0)
        y = solve_Ux_equals_y(u, z)
        if y is None:
            return None
        s = scalar_product(y, y)
        t = scalar_product(y, x0)
        x = divide_vector_by_scalar(y, vector_norm(y))
        lmbd = s / t
        lmbd_min = 1 / lmbd
        approx_pairs.append((lmbd_min, x))
        if abs(lmbd_min - lmbd0) <= eps:
            continue_iter = False
        else:
            x0 = x
            lmbd0 = lmbd_min
    return approx_pairs


def output_iter_approximations(pairs):
    if pairs is not None:
        print("Iteration\tEigenvalue\tEigenvector")
        for i, pair in enumerate(pairs):
            print(f"{i}: {round(pair[0], 5)}\t{pair[1]}")
    else:
        print("Not invertible")

