import numpy as np
import matplotlib.pyplot as plt

# ==========================
# Parâmetros do motor
# ==========================

Rs = 6.4
Rr = 3.5
Ls = 0.30
Lr = 0.30
Lm = 0.28

P = 4                 # número de polos
J = 0.02              # inércia
B = 0.001             # atrito viscoso
Tl = 0                # carga inicial

V_phase = 220/np.sqrt(3) * np.sqrt(2)
f = 60
we = 2*np.pi*f

# ==========================
# Simulação
# ==========================

dt = 1e-4
t_final = 2
N = int(t_final/dt)

t = np.linspace(0, t_final, N)

# Estados
psi_ds = 0
psi_qs = 0
psi_dr = 0
psi_qr = 0
wr = 0

# Vetores para salvar
wr_vec = []
Te_vec = []
ids_vec = []
ia_vec = []
ib_vec = []
ic_vec = []

den = Ls*Lr - Lm**2

for k in range(N):

    # Tensões trifásicas -> dq (partida direta)
    va = V_phase*np.sin(we*t[k])
    vb = V_phase*np.sin(we*t[k] - 2*np.pi/3)
    vc = V_phase*np.sin(we*t[k] + 2*np.pi/3)

    # Transformada de Clarke
    v_alpha = (2/3)*(va - 0.5*vb - 0.5*vc)
    v_beta  = (2/3)*(np.sqrt(3)/2*(vb - vc))

    # Park
    theta = we*t[k]
    vds =  v_alpha*np.cos(theta) + v_beta*np.sin(theta)
    vqs = -v_alpha*np.sin(theta) + v_beta*np.cos(theta)

    # Correntes a partir dos fluxos
    ids = (Lr*psi_ds - Lm*psi_dr)/den
    iqs = (Lr*psi_qs - Lm*psi_qr)/den
    idr = (-Lm*psi_ds + Ls*psi_dr)/den
    iqr = (-Lm*psi_qs + Ls*psi_qr)/den

    # dq -> alpha-beta
    i_alpha = ids*np.cos(theta) - iqs*np.sin(theta)
    i_beta  = ids*np.sin(theta) + iqs*np.cos(theta)

    # alpha-beta -> abc
    ia = i_alpha
    ib = -0.5*i_alpha + (np.sqrt(3)/2)*i_beta
    ic = -0.5*i_alpha - (np.sqrt(3)/2)*i_beta

    # Derivadas dos fluxos
    dpsi_ds = vds - Rs*ids + we*psi_qs
    dpsi_qs = vqs - Rs*iqs - we*psi_ds

    slip = we - (P/2)*wr

    dpsi_dr = -Rr*idr + slip*psi_qr
    dpsi_qr = -Rr*iqr - slip*psi_dr

    # Integração (Euler)
    psi_ds += dpsi_ds*dt
    psi_qs += dpsi_qs*dt
    psi_dr += dpsi_dr*dt
    psi_qr += dpsi_qr*dt

    # Torque
    Te = (3/2)*(P/2)*(psi_ds*iqs - psi_qs*ids)

    # Mecânica
    dwr = (Te - Tl - B*wr)/J
    wr += dwr*dt

    theta = we*t[k]   # ângulo elétrico síncrono

    # Salvar
    wr_vec.append(wr)
    Te_vec.append(Te)
    ids_vec.append(ids)
    ia_vec.append(ia)
    ib_vec.append(ib)
    ic_vec.append(ic)


# ==========================
# Gráficos
# ==========================

plt.figure()
plt.plot(t, wr_vec)
plt.title("Velocidade do Rotor (rad/s)")
plt.grid()

plt.figure()
plt.plot(t, Te_vec)
plt.title("Torque Eletromagnético")
plt.grid()

plt.figure()
plt.plot(t, ids_vec)
plt.title("Corrente Ids")
plt.grid()

plt.figure()
plt.plot(t, ia_vec, label="Ia")
plt.plot(t, ib_vec, label="Ib")
plt.plot(t, ic_vec, label="Ic")
plt.legend()
plt.grid()
plt.title("Correntes Trifásicas")

plt.show()
