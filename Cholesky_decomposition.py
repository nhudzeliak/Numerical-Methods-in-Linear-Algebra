import math


def get_U(A):
    n = len(A)
    U = [[None for i in range(n)] for j in range(n)]
    for i in range(n):
        for j in range(n):
            if j > i:
                temp_sum = 0
                for k in range(i):
                    temp_sum += U[k][i] * U[k][j]
                try:
                    U[i][j] = (A[i][j] - temp_sum) / U[i][i]
                except ZeroDivisionError:
                    return None
            elif j == i:
                temp_sum = 0
                for k in range(i):
                    temp_sum += U[k][i] ** 2
                try:
                    U[i][j] = math.sqrt(A[i][i] - temp_sum)
                except ValueError:
                    return None
            else:
                U[i][j] = 0
    return U


def solve(U, b):
    n = len(b)
    y = []
    for i in range(n):
        temp = 0
        for k in range(i):
            temp += U[k][i] * y[k]
        try:
            y.append((b[i]-temp) / U[i][i])
        except:
            return None
    x = [None for i in range(n)]
    for i in range(n-1, -1, -1):
        temp = 0
        for k in range(i+1, n):
            temp += U[i][k] * x[k]
        try:
            x[i] = (y[i] - temp) / U[i][i]
        except:
            return None
    x = [round(k, 4) for k in x]
    return x


def solution_to_string(x):
    if x is not None:
        s = ""
        for i in range(len(x)):
            s += f"x{i+1} = {x[i]}\n"
        return s
    else:
        return "Matrix A is not positive-definite."


def print_solution(x):
    print(solution_to_string(x))