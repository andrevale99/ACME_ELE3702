import numpy as np
from numpy import pi

# ============================================
# FUNCOES DE AUXILIO
# ============================================
def F_to_Ns(Poles, freqMotor):
    '''
    Funcao que retorna a velocidade
    de rotação do campo
    magnetomotriz

    Parameters:
        Poles: Quantidade de Polos do motor
        freqMotor: Frequencia do motor

    Returns:
        A conversão em RPM.
    '''
    return ( 2 / Poles * freqMotor * 60)

def Rad_to_Deg(rad):
    return (rad * 180. / np.pi)

def Deg_to_Rad(deg):
    return (deg * np.pi / 180.)
# ============================================
# 
# ============================================

class IM:

    # Vetor dos valores de escorregamento
    sVect = np.linspace(1e-3, 1, 1000)

    def __init__(self, ArrayParam=None, f=60, Vabc=None, 
                 _timeVect=np.empty(0), SimulationParam=[[0,0.1], 0.01]):
        if len(ArrayParam) < 13:
            print("Passar Parametros do Motor em um vetor de 8 posicoes\nna seguinte ordem:")
            print(f'Rs: Resistência do Estator\n \
                  Lms: Indutancia de magnetizacao do estator\n \
                  Lls: Indutancia de dispersao do Estator\n \
                  Ns: Voltas do Enrolamento do Estator\n \
                  Rr: Resitência do Rotor\n \
                  Lmr: Indutancia de magnetizacao do rotor\n \
                  Llr Indutancia do dispersao Rotor\n \
                  Nr: Voltas do Enrolamente do rotor\n \
                  Polos: Quantidade de Polos da maquina\n \
                  RPM: RPM maximo da maquina \n \
                  Xls: Reatancia do Estator \n \
                  Xm: Reatancia de magnetizacao \n \
                  Xlr: Reatancia do rotor \n \
                  ')
            print()
            print("Parametros passados independentes:")
            print(f'f: Frequencia da rede eletrica em Hz\n \
                  we: Frequencia da rede eletrica em rad/s\n \
                  Vabc: Vetor 1x3 Das tensoes das fases em RMS')
        else:
            self.Rs = ArrayParam[0]
            self.Lms = ArrayParam[1]
            self.Lls = ArrayParam[2]
            self.Ns = ArrayParam[3]

            self.Rr = ArrayParam[4]
            self.Lmr = ArrayParam[5]
            self.Llr = ArrayParam[6]
            self.Nr = ArrayParam[7]

            self.n = self.Nr / self.Ns

            self.Vabc = Vabc #Vetor 1x3 Das tensoes das fases em RMS
            self.f = f #Frequencia da Rede (Hz)
            self.we = 2 * pi * f #Frequencia da rede (rad/s)

            self.Polos = ArrayParam[8] # Polos da Maquina
            self.n = ArrayParam[9] #Velocidade Mecanica da maquina

            self.ws = self.f * 2 * pi * (2 / Polos) #Velocidade sincrona (rad/s)
            self.ns = 120 * self.f / Polos
            self.s = (self.ns - self.n) / self.ns

            self.Xls = ArrayParam[10]
            self.Xm = ArrayParam[11]
            self.Xlr = ArrayParam[12]

            # Calcula as reatancias caso sejam 0
            if self.Xls == 0:
                self.Xls = 2*pi*self.f*self.Lls
            if self.Xm == 0:
                self.Xm = 2*pi*self.f*self.Lms
            if self.Xlr == 0:
                self.Xlr = 2*pi*self.f*self.Llr

            # Impedancia do estator
            # Usando a simplificação pelo teorema
            # de Thevenin
            self.Zs = (1j * self.Xm * (self.Rs + 1j * self.Xls)) / (self.Rs + 1j * self.Xls + 1j*self.Xm)
            # Impedancia Rotor
            self.Zr = Rr / self.sVect + 1j*(self.Xlr)
            self.Z = self.Zs + self.Zr

            self.MtxRs = np.array([
                [self.Rs, 0, 0],
                [0, self.Rs, 0],
                [0, 0, self.Rs]
            ])

            self.MtxRr = np.array([
                [self.Rr, 0, 0],
                [0, self.Rr, 0],
                [0, 0, self.Rr]
            ])

            self.Ls = np.array([
                [(self.Lls + self.Lms), -self.Lms/2, -self.Lms/2],
                [-self.Lms/2, (self.Lls + self.Lms), -self.Lms/2],
                [-self.Lms/2, -self.Lms/2, (self.Lls + self.Lms)]
            ], dtype=np.float64)

            Lr = np.array([
                [(Lls + self.n**2 * Lms), -self.n**2 * Lms / 2, -self.n**2 * Lms / 2],
                [-self.n**2 * Lms / 2, (Lls + self.n**2 * Lms), -self.n**2 * Lms / 2],
                [-self.n**2 * Lms / 2, -self.n**2 * Lms / 2, (Lls + self.n**2 * Lms)]
            ], dtype=np.float64)


        self.t0 = SimulationParam[0][0]
        self.tf = SimulationParam[0][1]
        self.dt = SimulationParam[1]
        self.timeVect = 0
        # Trecho para criar um vetor de Tempo para simular
        # a posicao theta do rotor
        if _timeVect.size == 0:
            TemptimeVect = np.linspace(self.t0, self.tf, self.Vabc.shape[1])
            self.timeVect = np.copy(TemptimeVect) # Vetor de tempo temporario
            del TemptimeVect
        else:
            self.timeVect = np.copy(_timeVect)

    def IM_dqs(self, Fxyz=None, k=2/3, Theta=0):
        '''
        Funcao que retorna uma transformacao em quadratura 
        para os eixos dqn com o Referencial Estacionario (estator).

        Parametros:
        Fxyz: Vetor com as 3 variaveis (fluxo, Tensao, Corrente)
        k: Constante de magnitude (padrao=2/3)
        Theta: Angulo entre o eixo D com o enrolamento da fase A
        '''
        T =  k * np.array([
            [np.cos(Theta), np.cos(Theta - 2*pi/3), np.cos(Theta + 2*pi/3)],
            [-np.sin(Theta), -np.sin(Theta - 2*pi/3), -np.sin(Theta + 2*pi/3)],
            [1/2, 1/2, 1/2]
        ])

        return (T @ Fxyz)
    
    def IM_dqe(self, Fxyz=None, k=2/3):
        '''
        Funcao que retorna uma transformacao em quadratura 
        para os eixos dqn com o Referencial Rotatorio (rotor).

        Parametros:
        Fxyz: Vetor com as 3 variaveis (fluxo, Tensao, Corrente)
        k: Constante de magnitude (padrao=2/3)
        Theta: Angulo entre o eixo D com o enrolamento da fase A
        '''
        Vdqs = self.IM_dqs(Fxyz, k=k)

        Theta = self.we * self.timeVect

        R = np.array([
            [np.cos(Theta), np.sin(Theta)],
            [-np.sin(Theta), np.cos(Theta)],
        ])

        Vdqe = np.array([
            R[0,0,:] * Vdqs[0,:] + R[0,1,:] * Vdqs[1,:],
            R[1,0,:] * Vdqs[0,:] + R[1,1,:] * Vdqs[1,:],
            Vdqs[2,:]
        ])

        del R,Theta, Vdqs

        return Vdqe
    
    def IM_getIsIr(self):
        '''
        Retorna as correntes do estator (Is)
        e do rotor (Ir) do escorregamento maximo
        ate o minimo
        '''
        Is = np.max(self.Vabc[0]) / np.abs(self.Z)
        # Ismax = np.max(Is)

        # Corrente no rotor
        Ir = np.max(self.Vabc[0])  / (self.Rs + self.Rr / self.sVect + 1j*(self.Xls+self.Xlr))
        # Irmax =  np.max(self.Vabc[0])  / (self.Rs + self.Rr / 1 + 1j*(self.Xls+self.Xlr))

        return Is,Ir

    def IM_getTmec(self, Is, Ir):
        '''
        Retorna o torque mecanico do motor de indução
        considerando as 3 fases do motor e os valores
        escalres do Torque de Partida e Nominal.

        Parametros:
        Is: Corrente do estator
        Ir: Corrente do rotor
        '''


        Vs = np.max(self.Vabc[0])
        Tmec = (3 / self.ws) * ((Vs**2) / ((self.Rs + self.Rr / self.sVect) 
                                           + 1j*(self.Xls+self.Xlr))**2) * (self.Rr/self.sVect)
        
        sPartida = 1
        Tpartida = np.abs((3 / self.ws) * ((Vs**2) / ((self.Rs + self.Rr / sPartida) +
                                         1j*(self.Xls+self.Xlr))**2) * (self.Rr / sPartida))
        print(f'Tpartida = {Tpartida}\n')

        smax = self.Rr / np.sqrt(self.Rs**2 + (self.Xls+self.Xlr)**2)
        Tmax = np.abs((3 / self.ws) * ((Vs**2) / ((self.Rs + self.Rr / smax) +
                                     1j*(self.Xls+self.Xlr))**2) * (self.Rr / smax))

        Tnominal = (3 / self.ws) * ((Vs**2) / ((self.Rs + self.Rr / self.s) 
                                                         + 1j*(self.Xls+self.Xlr))**2) * (self.Rr/self.s)

        return Tmec, Tpartida, Tnominal


if __name__ == "__main__":
    import matplotlib.pyplot as plt

    #==========================================================
    # VARIAVEIS
    #==========================================================
    Samples = 1000
    timeVect = np.linspace(0, 0.15, Samples)

    Rs = 6.4 # Resistencia do estator
    Lls = 1.39e-3 # Indutancia de dispersao  do estator
    Lms = 41e-3 # Indutancia de magnetizacaodo estator
    Ns = 1 # Voltas no enrolamento do estator

    Rr = 4.25 # Resistencia do rotor
    Llr = 0.74e-3 # Indutancia de dispersao do rotor
    Lmr = 41e-3 # Indutancia de magnetizacao  do rotor
    Nr = 1 # Voltas no enrolamento do rotor

    Xls = 5.85 #Ohm
    Xlr = 5.85 #Ohm
    Xm = 137.8 #Ohm

    Vm = 220 # Amplitude da tensao da rede eletrica (V)

    f = 60 # Frequencia da rede (Hz)
    omega_e = 2 * pi * f # Frequencia da rede (rad/s)

    Polos = 4
    RPM = 1715

    Vabc = np.array([
        Vm*np.cos(omega_e*timeVect),
        Vm*np.cos(omega_e*timeVect - 2*np.pi/3),
        Vm*np.cos(omega_e*timeVect + 2*np.pi/3)
    ])

    #==========================================================

    Params = [Rs,Lms,Lls,Ns,Rr,Lmr,Llr,Nr,Polos,RPM,Xls,Xm,Xlr]
    x = IM(Params, Vabc=Vabc, f=f, _timeVect=timeVect)

    Vdqs = x.IM_dqs(Vabc)
    Vdqe = x.IM_dqe(Vabc)

    plt.title("Vdqn Park Transformation")
    plt.plot(timeVect, Vdqs[0,:], label='Vds')
    plt.plot(timeVect, Vdqs[1,:], label='Vqs')
    plt.plot(timeVect, Vdqe[0,:], label='Vde')
    plt.plot(timeVect, Vdqe[1,:], label='Vqe')
    plt.legend()
    plt.grid()
    plt.show()

    Is,Ir = x.IM_getIsIr()

    plt.title("Correntes")
    plt.plot((1-x.sVect)*x.ns, np.abs(Is), label="Is")
    plt.plot((1-x.sVect)*x.ns, np.abs(Ir), label="Ir")
    plt.legend()
    plt.grid()
    plt.show()

    Tmec = x.IM_getTmec(Is, Ir)
    plt.title("Torque Mecânico")
    plt.plot((1-x.sVect)*x.ns, np.abs(Tmec), label="Tmec")
    plt.legend()
    plt.grid()
    plt.show()