import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

#data1 = pd.read_csv("hindmarsh_rose_dopri5.csv")
data2 = pd.read_csv("hindmarsh_rose_dopri8.csv")
#data3 = pd.read_csv("hindmarsh_rose_rk4.csv")

fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d') 

#ax.plot(data1["Time"], data1["y"], data1["x"], label = "DOPRI5")
ax.plot(data2["Time"], data2["y"], data2["x"], color = 'red', label = 'DOPRI8')
#ax.plot(data3["Time"], data3["y"], data3["x"], color = 'green', label = 'RK4')

ax.set_xlabel("t")
ax.set_ylabel("y")
ax.set_zlabel("x")
ax.set_title("3D-график")

plt.show()