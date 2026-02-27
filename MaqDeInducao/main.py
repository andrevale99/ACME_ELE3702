from IMClass import IM, n_with_ns, tensao_vf
import numpy as np
from math import pi
import matplotlib.pyplot as plt


#==========================================================
# VARIAVEIS
#==========================================================
Samples = 1000
timeVect = np.linspace(0, 0.15, Samples)
f = 30 # Frequencia da rede (Hz)
omega_e = 2 * pi * f # Frequencia da rede (rad/s)

Vm = tensao_vf(f, f_min=10, f_max=60, V_min=50, V_max=220) # Amplitude da tensao da rede eletrica (V)

Polos = 4
RPM = 1715

Rs = 6.4 # Resistencia do estator
Lls = 1.39e-4 # Indutancia de dispersao  do estator
Lms = 41e-4 # Indutancia de magnetizacaodo estator
Ns = 1 # Voltas no enrolamento do estator
Xls = 5.85

Rr = 4.25 # Resistencia do rotor
Llr = 0.74e-4 # Indutancia de dispersao do rotor
Lmr = 41e-4 # Indutancia de magnetizacao  do rotor
Nr = 1 # Voltas no enrolamento do rotor
Xlr = 5.85

Xm = 137.8 #Ohm

Vabc = np.array([
    Vm*np.cos(omega_e*timeVect),
    Vm*np.cos(omega_e*timeVect - 2*np.pi/3),
    Vm*np.cos(omega_e*timeVect + 2*np.pi/3)
])
#==========================================================
Params = [Rs,Lms,Lls,Ns,Rr,Lmr,Llr,Nr,Polos,RPM,Xls,Xm,Xlr]
x = IM(Params, Vabc=Vabc, f=f, _timeVect=timeVect)

Tmec = x.IM_fsSpeedControl(50)
plt.plot(n_with_ns(x.sVect, x.ns), np.abs(Tmec), label="Tmec")
plt.legend()
plt.grid()
plt.show()

# Tmec1,_,_,_ = x1.IM_getTmec()
# Tmec2,_,_,_ = x2.IM_getTmec()
# plt.title("Torque Mecânico")
# plt.plot((1-x1.sVect)*x1.ns, np.abs(Tmec1), label="Tmec")
# plt.plot((1-x2.sVect)*x2.ns, np.abs(Tmec2), label="Tmec2")  
# plt.legend()
# plt.grid()
# plt.show()

# Vdqs = x.IM_dqs(Vabc)
# Vdqe = x.IM_dqe(Vabc)
# plt.title("Vdqn Park Transformation")
# plt.plot(timeVect, Vdqs[0,:], label='Vds')
# plt.plot(timeVect, Vdqs[1,:], label='Vqs')
# plt.plot(timeVect, Vdqe[0,:], label='Vde', ls='--')
# plt.plot(timeVect, Vdqe[1,:], label='Vqe', ls='--')
# plt.legend()
# plt.grid()
# plt.show()

# Is,Ir = x.IM_getIsIr()
# plt.title("Correntes")
# plt.plot((1-x.sVect)*x.ns, np.abs(Is), label="Is")
# plt.plot((1-x.sVect)*x.ns, np.abs(Ir), label="Ir")
# plt.legend()
# plt.grid()
# plt.show()

# Tmec,Tp,Tn,[Tmax,smax] = x.IM_getTmec()
# plt.title("Torque Mecânico")
# plt.plot((1-x.sVect)*x.ns, np.abs(Tmec), label="Tmec")
# plt.scatter((1-smax)*x.ns, Tmax, label="Tmec")
# plt.legend()
# plt.grid()
# plt.show()

# print(x.IM_getIl(15, 208, 'F'))