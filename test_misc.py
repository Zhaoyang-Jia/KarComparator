from scipy.optimize import linear_sum_assignment
import numpy as np

cost = np.array([[1, 2, 3], [1, 2, 3]])
row_ind, col_ind = linear_sum_assignment(cost)
print(row_ind)
print(col_ind)
