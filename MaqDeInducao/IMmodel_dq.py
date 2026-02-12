import matplotlib.pyplot as plt
import numpy as np
from numpy import pi

dt = 0.001

t = np.arange(0, 0.8, dt)

Vm = 220
f = 100 # Frequencia da rede (Hz)
omega_e = 2 * pi * f # Frequencia da rede (rad/s)

Vabc = np.array([
    Vm*np.cos(omega_e*t),
    Vm*np.cos(omega_e*t - 2*np.pi/3),
    Vm*np.cos(omega_e*t + 2*np.pi/3)
])

Theta = 0
T =  2/3 * np.array([
    [np.cos(Theta), np.cos(Theta - 2*pi/3), np.cos(Theta + 2*pi/3)],
    [-np.sin(Theta), -np.sin(Theta - 2*pi/3), -np.sin(Theta + 2*pi/3)],
    [1/2, 1/2, 1/2]
])

Vdqs = T @ Vabc

Rs = 1

lambda_s_ds = np.zeros((Vdqs.shape[1]))
for k in range(N):

    derivada = vds[k] - Rs * ids[k]
    lambda_s_ds[k] = lambda_s_ds[k] + derivada * dt