import numpy as np
from numpy import pi

class IM:

    def __init__(self, ArrayParam=None, f=60,Vabc=None):
        if ArrayParam == None:
            print("Passar Parametros do Motor em um vetor de 8 posicoes\nna seguinte ordem:")
            print(f'Rs: Resistência do Estator\n \
                  Lms: Indutancia de magnetizacao do estator\n \
                  Lls: Indutancia de dispersao do Estator\n \
                  Ns: Voltas do Enrolamento do Estator\n \
                  Rr: Resitência do Rotor\n \
                  Lmr: Indutancia de magnetizacao do rotor\n \
                  Llr Indutancia do dispersao Rotor\n \
                  Nr: Voltas do Enrolamente do rotor\n \
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

            self.Rs = np.array([
                [self.Rs, 0, 0],
                [0, self.Rs, 0],
                [0, 0, self.Rs]
            ])


            self.Rr = np.array([
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

    def IM_dqs(self, Fxyz=None ,k=2/3,Theta=0):
        '''
        Funcao que retorna uma transformacao em quadratura 
        para os eixos dqn.

        Parametros:
        Fxyz: Vetor com as 3 variaveis
        k: Constante de magnitude (padrao=2/3)
        Theta: Angulo entre o eixo D com o enrolamento da fase A
        '''
        T =  k * np.array([
            [np.cos(Theta), np.cos(Theta - 2*pi/3), np.cos(Theta + 2*pi/3)],
            [-np.sin(Theta), -np.sin(Theta - 2*pi/3), -np.sin(Theta + 2*pi/3)],
            [1/2, 1/2, 1/2]
        ])

        return (T @ Fxyz)
    
    def IM_dqe(self, Fxyz=None , k=2/3, timeVect=np.empty(0)):

        if timeVect.size == 0:
            _timeVect = np.linspace(0, 0.015, Fxyz.shape[1])
            timeVect = np.copy(_timeVect) # Vetor de tempo temporario
            del _timeVect

        Vdqs = self.IM_dqs(Fxyz)

        Theta = self.we * timeVect

        R = np.array([
            [np.cos(Theta), np.sin(Theta)],
            [-np.sin(Theta), np.cos(Theta)]
        ])

        Vdqe = np.array([
            R[0,0,:] * Vdqs[0,:] + R[0,1,:] * Vdqs[1,:],
            R[1,0,:] * Vdqs[0,:] + R[1,1,:] * Vdqs[1,:]
        ])


        del R,Theta, Vdqs

        return Vdqe


if __name__ == "__main__":
    import matplotlib.pyplot as plt

    #==========================================================
    # VARIAVEIS
    #==========================================================
    Samples = 1000
    timeVect = np.linspace(0, 0.015, Samples)

    Rs = 0.294 # Resistencia do estator
    Lls = 1.39e-3 # Indutancia de dispersao  do estator
    Lms = 41e-3 # Indutancia de magnetizacaodo estator
    Ns = 1 # Voltas no enrolamento do estator

    Rr = 0.156 # Resistencia do rotor
    Llr = 0.74e-3 # Indutancia de dispersao do rotor
    Lmr = 41e-3 # Indutancia de magnetizacao  do rotor
    Nr = 1 # Voltas no enrolamento do rotor

    Vm = 100 # Amplitude da tensao da rede eletrica (V)

    f = 100 # Frequencia da rede (Hz)
    omega_e = 2 * pi * f # Frequencia da rede (rad/s)

    Vabc = np.array([
        Vm*np.cos(omega_e*timeVect),
        Vm*np.cos(omega_e*timeVect - 2*np.pi/3),
        Vm*np.cos(omega_e*timeVect + 2*np.pi/3)
    ])

    #==========================================================

    Params = [Rs,Lms,Lls,Ns,Rr,Lmr,Llr,Nr]
    x = IM(Params,Vabc=Vabc, f=f)

    Vdqs = x.IM_dqs(Vabc)

    plt.title("Vdqn estacionario (referencia no estator)")
    plt.plot(timeVect, Vdqs[0,:], label='Vds')
    plt.plot(timeVect, Vdqs[1,:], label='Vqs')
    plt.legend()
    plt.show()

    Vdqe = x.IM_dqe(Vabc, timeVect=timeVect)

    plt.title("Vdqn rotacional (referencia no rotor)")
    plt.plot(timeVect, Vdqe[0,:], label='Vde')
    plt.plot(timeVect, Vdqe[1,:], label='Vqe')
    plt.legend()
    plt.show()