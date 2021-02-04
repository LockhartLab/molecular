import matplotlib.pyplot as plt
import molecular as mol

n_temps = 10
geometric = mol.temp_schedule(300, 440, n_temps, 'geometric')
linear = mol.temp_schedule(300, 440, n_temps, 'linear')
parabolic = mol.temp_schedule(300, 440, n_temps, 'parabolic')

plt.figure()
plt.plot(range(n_temps), geometric, label='geometric')
plt.plot(range(n_temps), linear, label='linear')
plt.plot(range(n_temps), parabolic, label='parabolic')
plt.xlabel('index')
plt.ylabel('temperature')
plt.legend()
plt.show()