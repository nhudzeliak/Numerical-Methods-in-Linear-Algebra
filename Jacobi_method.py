def check_diagonal_domination(A):
    dominant = True
    for i in range(len(A)):
        sum_row_no_diag = 0
        for j in range(len(A)):
            if i != j:
                sum_row_no_diag += abs(A[i][j])
        if abs(A[i][i]) <= sum_row_no_diag:
            dominant = False
            return dominant
    return dominant


def normalize(A, b):
    c = []
    d = []
    for i in range(len(A)):
        c_row = []
        for j in range(len(A)):
            if i == j:
                c_row.append(0)
            else:
                c_row.append(-A[i][j] / A[i][i])
        c.append(c_row)
        d.append(b[i] / A[i][i])
    return c, d


def norm_matr(c):
    sums = []
    for row in c:
        sum_abs_row = 0
        for j in range(len(c)):
            sum_abs_row += abs(row[j])
        sums.append(sum_abs_row)
    norm_c = max(sums)
    return norm_c


def Jacobi(A, b, x_prev, e):
    c, d = normalize(A, b)
    if check_diagonal_domination(A):
        norm_c = norm_matr(c)
        continue_iter = True
        all_values = [x_prev]
        while continue_iter:
            x = []
            for i in range(len(A)):
                temp_sum = 0
                for j in range(len(A)):
                    if i != j:
                        temp_sum += A[i][j] * x_prev[j] / A[i][i]
                x.append((-1) * temp_sum + b[i] / A[i][i])
            all_values.append(x)
            subtr_x = [abs(x[i] - x_prev[i]) for i in range(len(x))]
            if norm_c * max(subtr_x) / (1 - norm_c) < e:
                continue_iter = False
            x_prev = x
        return all_values
    else:
        print("Not diagonally dominant")


def iterations_output(iter_values):
    for i, row in enumerate(iter_values):
        print(f"{i}: {row}")


