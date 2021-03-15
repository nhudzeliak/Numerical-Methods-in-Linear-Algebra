def Inverse(A, n):
    """Returns inverse matrix implementing Gaussian elimination.
    :param A: input square matrix to inverse
    :param n: size of the matrix A
    :return: inverse of matrix A
    """
    x = [[None for i in range(n)] for j in range(n)]
    b = [[0 for i in range(n)] for j in range(n)]
    for p in range(n):
        b[p][p] = 1

    eps = 0.000001
    for k in range(0, n-1):
        max_m = k
        for m in range(k+1, n):
            if abs(A[m][k]) > abs(A[max_m][k]):
                max_m = m
        if abs(A[max_m][k]) < eps:
            return None
        else:
            A[max_m], A[k] = A[k], A[max_m]
            b[max_m], b[k] = b[k], b[max_m]

        for i in range(k+1, n):
            m = (-1) * A[i][k] / A[k][k]
            for j in range(k+1, n):
                A[i][j] = A[i][j] + m * A[k][j]
            for j in range(0, n):
                b[i][j] = b[i][j] + m * b[k][j]

    if abs(A[n-1][n-1]) < eps:
        return None

    for j in range(0, n):
        x[n-1][j] = (b[n-1][j] / A[n-1][n-1])
    for k in range(n-2, -1, -1):
        for i in range(n):
            temp = 0
            for j in range(k+1, n):
                temp += A[k][j] * x[j][i]
            x[k][i] = (b[k][i]-temp) / A[k][k]
    return x

