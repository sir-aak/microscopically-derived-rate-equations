import numpy as np
# ~ import matplotlib.pyplot as plt


A = np.array([[0.0, 1.0, 1.0, 1.0], 
              [1.0, 0.0, 1.0, 1.0], 
              [1.0, 1.0, 0.0, 1.0], 
              [1.0, 1.0, 1.0, 0.0]])

# ~ l1 = 3.0
# ~ l2 = -1.0
# ~ l3 = -1.0
# ~ l4 = -1.0

# ~ EW = np.array([l1, l2, l3, l4])

ew, S = np.linalg.eig(A)

D = np.diag(ew)

print(S @ D @ np.linalg.inv(S))
