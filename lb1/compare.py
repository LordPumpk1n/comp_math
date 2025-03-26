import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

data_dopri5 = pd.read_csv("hindmarsh_rose_dopri5.csv")
data_dopri8 = pd.read_csv("hindmarsh_rose_dopri8.csv")
data_rk4 = pd.read_csv("hindmarsh_rose_rk4.csv")

time = data_dopri8["Time"]

error_dopri5_x = np.abs(data_dopri5["x"] - data_dopri8["x"])
error_dopri5_y = np.abs(data_dopri5["y"] - data_dopri8["x"])

error_rk4_x = np.abs(data_rk4["x"] - data_dopri8["x"])
error_rk4_y = np.abs(data_rk4["y"] - data_dopri8["x"])

plt.figure(figsize=(12, 6))

plt.subplot(1, 2, 1)
plt.plot(time, error_dopri5_x, label="DOPRI5 vs DOPRI8 (x)", color="blue")
plt.plot(time, error_rk4_x, label="RK4 vs DOPRI8 (x)", color="green")
plt.xlabel("Time (t)")
plt.ylabel("Absolute Error (x)")
plt.yscale("log")  # Логарифмическая шкала для наглядности
plt.title("Ошибка по переменной x")
plt.legend()

plt.subplot(1, 2, 2)
plt.plot(time, error_dopri5_y, label="DOPRI5 vs DOPRI8 (y)", color="red")
plt.plot(time, error_rk4_y, label="RK4 vs DOPRI8 (y)", color="orange")
plt.xlabel("Time (t)")
plt.ylabel("Absolute Error (y)")
plt.yscale("log")
plt.title("Ошибка по переменной y")
plt.legend()

plt.tight_layout()
plt.show()