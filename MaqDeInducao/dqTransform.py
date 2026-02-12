import matplotlib.pyplot as plt
import numpy as np
from numpy import pi

Samples = 1000

timeVect = np.linspace(0, 0.015, Samples)

Vm = 100 # Amplitude da tensao da rede eletrica (V)

f = 100 # Frequencia da rede (Hz)
omega_e = 2 * pi * f # Frequencia da rede (rad/s)

Vabc = np.array([
    Vm*np.cos(omega_e*timeVect),
    Vm*np.cos(omega_e*timeVect - 2*np.pi/3),
    Vm*np.cos(omega_e*timeVect + 2*np.pi/3)
])

Theta = 0
T =  2/3 * np.array([
    [np.cos(Theta), np.cos(Theta - 2*pi/3), np.cos(Theta + 2*pi/3)],
    [-np.sin(Theta), -np.sin(Theta - 2*pi/3), -np.sin(Theta + 2*pi/3)],
    [1/2, 1/2, 1/2]
])


Vdqs = T @ Vabc

Theta = omega_e*timeVect

R = np.array([
    [np.cos(Theta), np.sin(Theta)],
    [-np.sin(Theta), np.cos(Theta)]
])

Vdqe = np.array([
    R[0,0,:] * Vdqs[0,:] + R[0,1,:] * Vdqs[1,:],
    R[1,0,:] * Vdqs[0,:] + R[1,1,:] * Vdqs[1,:]
])

plt.subplot(211)

plt.title("Vabc estator")
plt.plot(timeVect, Vabc[0,:], label='Va')
plt.plot(timeVect, Vabc[1,:], label='Vb')
plt.plot(timeVect, Vabc[2,:], label='Vc')
plt.grid()
plt.legend()

plt.subplot(212)
plt.title("Vdqn estacionario (referencia no estator)")
plt.plot(timeVect, Vdqs[0,:], label='Vds')
plt.plot(timeVect, Vdqs[1,:], label='Vqs')
plt.plot(timeVect, Vdqe[0,:], label='Vde', linestyle='--')
plt.plot(timeVect, Vdqe[1,:], label='Vqe', linestyle='-.')
plt.grid()
plt.legend()
plt.tight_layout()

plt.show()