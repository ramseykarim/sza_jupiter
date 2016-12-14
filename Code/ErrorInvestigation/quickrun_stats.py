import unpack as up
import numpy as np
import matplotlib.pyplot as plt

u = up.Unpack()
u.error_investigation.execute_testing(100)
res = u.error_investigation.results[0].results
xs = np.arange(res[0].size) + 2.

plt.figure(1)
plt.title("Error Averages")
n_plots = 5
for i in range(n_plots):
    plt.plot(xs, res[i + 2])
plt.xlabel("Number of points used")
plt.ylabel("Error (sort of jy)")
plt.legend([str(x + 1) for x in range(n_plots)])
plt.yscale('log')
plt.show()
