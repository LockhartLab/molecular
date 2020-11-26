
from scipy.constants import Boltzmann, Avogadro, atmosphere

# Boltzmann constant in kcal/K/mol
kb = Boltzmann * Avogadro / 4184.

# Pressure
pressure = Avogadro * atmosphere / 1e30 / 4184.

if __name__ == '__main__':
    print(f'kb={kb} kcal/mol')