import matplotlib.pyplot as plt
import numpy as np
from numpy import pi

Samples = 1000
timeVect = np.linspace(0, 0.015, Samples)

Vm = 100
f = 100
omega_e = 2 * pi * f

# -------------------------
# Tensões trifásicas
# -------------------------
Vabc = np.array([
    Vm*np.cos(omega_e*timeVect),
    Vm*np.cos(omega_e*timeVect - 2*pi/3),
    Vm*np.cos(omega_e*timeVect + 2*pi/3)
])

# -------------------------
# Clarke (abc → αβ0)
# -------------------------
T_clarke = (2/3) * np.array([
    [1, -1/2, -1/2],
    [0, np.sqrt(3)/2, -np.sqrt(3)/2],
    [1/2, 1/2, 1/2]
])

Valpha_beta_0 = T_clarke @ Vabc

Valpha = Valpha_beta_0[0,:]
Vbeta  = Valpha_beta_0[1,:]
print(Vabc.shape)
print(Valpha.shape)

# -------------------------
# Ângulo do rotor
# -------------------------
omega_r = omega_e  # rotor mais lento que campo
theta = omega_r * timeVect
print(theta.shape)

# -------------------------
# Park (αβ → dq rotor)
# -------------------------
Vd =  Valpha*np.cos(theta) + Vbeta*np.sin(theta)
Vq = -Valpha*np.sin(theta) + Vbeta*np.cos(theta)

# -------------------------
# Plot
# -------------------------
plt.subplot(211)
plt.title("Vabc - Estator")
plt.plot(timeVect, Vabc[0,:], label='Va')
plt.plot(timeVect, Vabc[1,:], label='Vb')
plt.plot(timeVect, Vabc[2,:], label='Vc')
plt.grid()
plt.legend()

plt.subplot(212)
plt.title("Transformações")
plt.plot(timeVect, Valpha, label='Valpha')
plt.plot(timeVect, Vbeta, label='Vbeta')
plt.plot(timeVect, Vd, '--', label='Vd (rotor)')
plt.plot(timeVect, Vq, '-.', label='Vq (rotor)')
plt.grid()
plt.legend()

plt.tight_layout()
plt.show()
