import unpack as up
import numpy as np
import matplotlib.pyplot as plt

u = up.Unpack()
u.error_investigation.execute_testing(100)
res = u.error_investigation.results[0].results
xs = np.arange(res[0].size) + 2.

plt.figure(1)
plt.title("Std Dev")
plt.plot(xs, res[0])
plt.xlabel("Number of points used")
plt.ylabel("Std Dev (jy)")

plt.figure(3)
plt.title("Error Averages")
for i in range(4):
    plt.plot(xs, res[i + 2])
plt.xlabel("Number of points used")
plt.ylabel("Error (sort of jy)")
plt.legend(['1', '2', '3', '4'])
plt.yscale('log')
plt.show()

"""
plt.figure(0)
plt.xlabel("Frequency (GHz)")
plt.ylabel("$T_B$ (K)")
plt.title("One of the average formulas")
u.plot_it()
f = open("formulation1.txt", 'w')
for ch in u.channel_obj_list:
    f.write(str(ch.tb) + ', ' + str(ch.tb_error) + '\n')
f.close()
"""
