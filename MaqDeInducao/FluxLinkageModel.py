import numpy as np
import matplotlib.pyplot as plt

'''
Modelagem de um motor de inducao a partir
do fluxo magnetico dos enrolamentos do motor
'''

# Variaveis Globais do sistema

rs = 10 # Resistencia do estator
Lls = 1 # Indutancia de magnetizacao do estator
Lms = 1 # Indutancia de dispersao do estator
Ns = 1 # Voltas no enrolamento do estator

Rs = np.array([
    [rs, 0, 0],
    [0, rs, 0],
    [0, 0, rs]
])

rr = 1 # Resistencia do rotor
Llr = 1 # Indutancia de magnetizacao do rotor
Lmr = 1 # Indutancia de dispersao do rotor
Nr = 2 # Voltas no enrolamento do rotor

Rr = np.array([
    [rr, 0, 0],
    [0, rr, 0],
    [0, 0, rr]
])

n = Nr / Ns # Razao de espiras

Samples = 1000# Quantidade de dados para a simulacao

# Posicao do rotor
theta_r = np.linspace(0, np.pi, Samples)
# theta_r = 0

t = np.linspace(0,0.05,Samples) # Vetor tempo para simulacao
# t = 0 

f = 60 # Frequencia da rede eletrica em Hz

omega = 2 * np.pi * f # Frequencia da rede em rad/s

I = 1 # Amplitude da correte no estator

# Correntes nas fases
Iabcs = np.array([
    I*np.cos(omega*t),
    I*np.cos(omega*t - 2*np.pi/3),
    I*np.cos(omega*t + 2*np.pi/3)
], dtype=np.float64)

# Vetor temporario para realizar os calculos
Iabcr = np.ones((3, Samples)) * (I/2)

# Diagonal principal = indutancia do estator do enrolamento
# da fase. Sendo a soma das indutancias de magnetizacao e de
# dispersão.
# Triangular Inferior/Superior = indutancia de magnetizacao 
# entre as fases dos enrolamentos. Como estão defasadas em 120 graus
# Lasbs = Lascs = Lbsas = Lbscs = Lcsas = Lcsbs = cos(2pi/3)Lms

# Considerando a distancia entreferro constate, as indutancias 
# possuem valores constantes
Ls = np.array([
    [(Lls + Lms), -Lms/2, -Lms/2],
    [-Lms/2, (Lls + Lms), -Lms/2],
    [-Lms/2, -Lms/2, (Lls + Lms)]
], dtype=np.float64)

# A logica da matriz de indutancia do rotor é a mesma
# da matriz de indutancias do estator, porem, na do rotor
# é considerado a razao entre o numero de voltas do ensolamento
# do rotor e estator nas indutancias de magnetizacao da fase e 
# entre as fases.
# Larbr = Larcr = Lbrar = L brcr = Lcrbr = Lmr cos(2pi/3) = -1/2 Lmr
#  = -1/2 (Nr/Ns)^2 Lms
Lr = np.array([
    [(Lls + n**2 * Lms), -n**2 * Lms / 2, -n**2 * Lms / 2],
    [-n**2 * Lms / 2, (Lls + n**2 * Lms), -n**2 * Lms / 2],
    [-n**2 * Lms / 2, -n**2 * Lms / 2, (Lls + n**2 * Lms)]
], dtype=np.float64)

# A matriz de induntacia entre as fases do estator e do rotor 
# eh modelada a partir da posicao do rotor com o eixo definido
# a partir de uma referencia do estator theta_r, e da razao de 
# espiras de enrolamento.
# Lxsyr = Nr/Ns Lms cos(theta_r)
# sendo x = a,b,c e y=a,b,c
# Lsr = np.array([
#     [np.cos(theta_r), np.cos(theta_r + 2*np.pi/3), np.cos(theta_r - 2*np.pi/3)],
#     [np.cos(theta_r - 2*np.pi/3), np.cos(theta_r),np.cos(theta_r + 2*np.pi/3)],
#     [np.cos(theta_r + 2*np.pi/3), np.cos(theta_r - 2*np.pi/3), np.cos(theta_r)]
# ], dtype=np.float64)

lambda_s = np.zeros((3, Samples))
lambda_r = np.zeros((3, Samples))
for k in range(Samples):

    # Matriz Lsr no instante k
    Lsr = np.array([
        [np.cos(theta_r[k]), np.cos(theta_r[k] + 2*np.pi/3), np.cos(theta_r[k] - 2*np.pi/3)],
        [np.cos(theta_r[k] - 2*np.pi/3), np.cos(theta_r[k]), np.cos(theta_r[k] + 2*np.pi/3)],
        [np.cos(theta_r[k] + 2*np.pi/3), np.cos(theta_r[k] - 2*np.pi/3), np.cos(theta_r[k])]
    ])

    # Fluxo estator
    lambda_s[:, k] = (
        Ls @ Iabcs[:, k]
        + n * Lms * Lsr @ Iabcr[:, k]
    )

    # Fluxo rotor
    lambda_r[:, k] = (
        Lr @ Iabcr[:, k]
        + n * Lms * Lsr.T @ Iabcs[:, k]
    )

# O fluxo no estator é a soma das induntacias do estator
# com as correntes do estator mais a indutancias entre o estator e 
# rotor com as correntes do rotor
# lambda_abcs = Ls*Iabcs + Lsr*Iabcr = Ls*Iabcs + n*Lms*Lsr*Iabcr
# 
# O mesmo valor para o rotor
# lambda_abcr = Lr*Iabcr + Lrs*Iabcs = Lr*Iabcr + n*Lms*Lsr^T*Iabcs

# O fluxo total do motor de inducao (fluxos do estator e rotor),
# pode ser calculado por:
# lambda_as = lambda_asas * lambda_asbs + 
#             lambda_ascs + lambda_asbr + lambda_ascr 
# ou
# lambda_as = Lasas*Ias + Lasbs*Ibs + Lascs*Ics + 
#               Lasar*Iar + Lasbr*Ibr + Lascr*Icr
# 
# A tensao das fases Vabcs no estator e calculada
# substituindo os fluxos concatenados pela equacao das
# repectivas tensoes em suas fases
# Vabcs = Rs * Iabcs + d/dt(lambda_abcs)
#       =  Rs * Iabcs + Ls * d/dt(Iabcs) + Lsr * d/dt(Iabcs)


d_lambda_s = np.gradient(lambda_s, t, axis=1)
Vabcs = Rs @ Iabcs + d_lambda_s

plt.plot(t, Vabcs[0,:], label="Va")
plt.plot(t, Vabcs[1,:], label="Vb")
plt.plot(t, Vabcs[2,:], label="Vc")
plt.grid()
plt.legend()
plt.show()