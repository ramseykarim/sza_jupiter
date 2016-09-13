import unpack as up
import matplotlib.pyplot as plt

u = up.Unpack()
plt.figure(0)
plt.xlabel("Frequency (GHz)")
plt.ylabel("$T_B$ (K)")
plt.title("One of the average formulas")
u.plotit()
f = open("formulation1.txt", 'w')
for ch in u.channel_obj_list:
    f.write(str(ch.tb) + ', ' + str(ch.tb_error) + '\n')
f.close()

