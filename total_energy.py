import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import CubicSpline


g = 9.81
c = 2/5
m = 0.031

h = 0.200
xfast=np.asarray([0,h,2*h,3*h,4*h,5*h,6*h,7*h])

x_values_1 = np.genfromtxt("./data/1_data.txt", delimiter=";", skip_header=2, usecols=1)
x_values_2 = np.genfromtxt("./data/2_data.txt", delimiter=";", skip_header=2, usecols=1)
x_values_3 = np.genfromtxt("./data/3_data.txt", delimiter=";", skip_header=2, usecols=1)
x_values_4 = np.genfromtxt("./data/4_data.txt", delimiter=";", skip_header=2, usecols=1)
x_values_5 = np.genfromtxt("./data/5_data.txt", delimiter=";", skip_header=2, usecols=1)
x_values_6 = np.genfromtxt("./data/6_data.txt", delimiter=";", skip_header=2, usecols=1)
x_values_7 = np.genfromtxt("./data/7_data.txt", delimiter=";", skip_header=2, usecols=1)
x_values_8 = np.genfromtxt("./data/8_data.txt", delimiter=";", skip_header=2, usecols=1)
x_values_9 = np.genfromtxt("./data/9_data.txt", delimiter=";", skip_header=2, usecols=1)
x_values_10 = np.genfromtxt("./data/10_data.txt", delimiter=";", skip_header=2, usecols=1)


y_values1 = np.genfromtxt("./data/1_data.txt", delimiter=";", skip_header=2, usecols=2)
v1 = np.genfromtxt("./data/1_data.txt", delimiter=";", skip_header=2, usecols=5)
y_values2 = np.genfromtxt("./data/2_data.txt", delimiter=";", skip_header=2, usecols=2)
v2 = np.genfromtxt("./data/2_data.txt", delimiter=";", skip_header=2, usecols=5)
y_values3 = np.genfromtxt("./data/3_data.txt", delimiter=";", skip_header=2, usecols=2)
v3 = np.genfromtxt("./data/3_data.txt", delimiter=";", skip_header=2, usecols=5)
y_values4 = np.genfromtxt("./data/4_data.txt", delimiter=";", skip_header=2, usecols=2)
v4 = np.genfromtxt("./data/4_data.txt", delimiter=";", skip_header=2, usecols=5)
y_values5 = np.genfromtxt("./data/5_data.txt", delimiter=";", skip_header=2, usecols=2)
v5 = np.genfromtxt("./data/5_data.txt", delimiter=";", skip_header=2, usecols=5)
y_values6 = np.genfromtxt("./data/6_data.txt", delimiter=";", skip_header=2, usecols=2)
v6 = np.genfromtxt("./data/6_data.txt", delimiter=";", skip_header=2, usecols=5)
y_values7 = np.genfromtxt("./data/7_data.txt", delimiter=";", skip_header=2, usecols=2)
v7 = np.genfromtxt("./data/7_data.txt", delimiter=";", skip_header=2, usecols=5)
y_values8 = np.genfromtxt("./data/8_data.txt", delimiter=";", skip_header=2, usecols=2)
v8 = np.genfromtxt("./data/8_data.txt", delimiter=";", skip_header=2, usecols=5)
y_values9 = np.genfromtxt("./data/9_data.txt", delimiter=";", skip_header=2, usecols=2)
v9 = np.genfromtxt("./data/9_data.txt", delimiter=";", skip_header=2, usecols=5)
y_values10 = np.genfromtxt("./data/10_data.txt", delimiter=";", skip_header=2, usecols=2)
v10 = np.genfromtxt("./data/10_data.txt", delimiter=";", skip_header=2, usecols=5)

def total(y_values, v):
  h = []
  for value in y_values:
    h.append(value - min(y_values))

  v[0] = 0

  #Her lager vi en tabell med x-verdier mellom 0 og 1.4 m
  xmin = 0.000
  xmax = 1.401
  dx = 0.001
  x = np.arange(xmin, xmax, dx)   

  def total_energy(v, h):
    new_list = []

    for i in range(len(v)):
      # print("kinetic: ", (1+c)*0.5*m*v[i]**2)
      # print("height: ", m*g*h[i])
      energy = (1+c)*0.5*m*(v[i]**2) + m*g*h[i]
      new_list.append(energy)

    return new_list

  total_energy_list = total_energy(v, h)
  total_energy_start = m*g*h[0]
  return total_energy_list, total_energy_start


total_energy_list_1, total_energy_start_1 = total(y_values1, v1)
total_energy_list_2, total_energy_start_2 = total(y_values2, v2)
total_energy_list_3, total_energy_start_3 = total(y_values3, v3)
total_energy_list_4, total_energy_start_4 = total(y_values4, v4)
total_energy_list_5, total_energy_start_5 = total(y_values5, v5)
total_energy_list_6, total_energy_start_6 = total(y_values6, v6)
total_energy_list_7, total_energy_start_7 = total(y_values7, v7)
total_energy_list_8, total_energy_start_8 = total(y_values8, v8)
total_energy_list_9, total_energy_start_9 = total(y_values9, v9)
total_energy_list_10, total_energy_start_10 = total(y_values10, v10)


baneform = plt.figure('Graf',figsize=(12,6))

# plt.plot(x_values_1, [total_energy_start_1]*len(x_values_1))
# plt.plot(x_values_1, total_energy_list_1)
# plt.plot(x_values_2, [total_energy_start_2]*len(x_values_2))
# plt.plot(x_values_2, total_energy_list_2)
plt.plot(x_values_3, [total_energy_start_3]*len(x_values_3))
plt.plot(x_values_3, total_energy_list_3)
# plt.plot(x_values_4, [total_energy_start_4]*len(x_values_4))
# plt.plot(x_values_4, total_energy_list_4)
plt.plot(x_values_5, [total_energy_start_5]*len(x_values_5))
plt.plot(x_values_5, total_energy_list_5)
# plt.plot(x_values_6, [total_energy_start_6]*len(x_values_6))
# plt.plot(x_values_6, total_energy_list_6)
# plt.plot(x_values_7, [total_energy_start_7]*len(x_values_7))
# plt.plot(x_values_7, total_energy_list_7)
# plt.plot(x_values_8, [total_energy_start_8]*len(x_values_8))
# plt.plot(x_values_8, total_energy_list_8)
# plt.plot(x_values_9, [total_energy_start_9]*len(x_values_9))
# plt.plot(x_values_9, total_energy_list_9)
# plt.plot(x_values_10, [total_energy_start_10]*len(x_values_10))
# plt.plot(x_values_10, total_energy_list_10)

plt.title('Total energi langs banen')
plt.xlabel('$x$ (m)',fontsize=20)
plt.ylabel('$E_T$ (J)',fontsize=20)
# plt.ylim(-1,2.5)
plt.grid()
plt.show()