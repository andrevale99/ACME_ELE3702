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

def n_with_ns(sVect, ns):
    '''
    Retorna a velocidade mecanica a partir da
    velocidade sincrona e do escorregamenta

    Parametro:
    sVect: Vetor de escorregamento [0,1]
    ns: Velocidade sincrona
    '''
    return (1-sVect)*ns

def tensao_vf(f, f_min=10, f_max=60, V_min=50, V_max=220):
    """
    Calcula tensão linear em função da frequência (controle V/f)

    f      : frequência (Hz) (escalar ou array)
    f_min  : frequência mínima
    f_max  : frequência máxima
    V_min  : tensão mínima
    V_max  : tensão máxima
    """

    # Limita a frequência dentro do intervalo
    f = np.clip(f, f_min, f_max)

    # Lei linear
    V = V_min + (f - f_min) * (V_max - V_min) / (f_max - f_min)

    return V
# ============================================
# 
# ============================================

class IM:

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

            self.Nsr = self.Nr / self.Ns # Razao entre os enrolamentos do rotor e do estator

            self.Vabc = Vabc #Vetor 1x3 Das tensoes das fases em RMS
            self.f = f #Frequencia da Rede (Hz)
            self.we = 2 * pi * f #Frequencia da rede (rad/s)

            self.Polos = ArrayParam[8] # Polos da Maquina
            self.n = ArrayParam[9] #Velocidade Mecanica da maquina

            self.ws = self.f * 2 * pi * (2 / self.Polos) #Velocidade sincrona (rad/s)
            self.ns = 120 * self.f / self. Polos
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
            self.Zr = self.Rr / self.sVect + 1j*(self.Xlr)
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
                [(self.Lls + self.Nsr**2 * self.Lms), -self.Nsr**2 * self.Lms / 2, -self.Nsr**2 * self.Lms / 2],
                [-self.Nsr * self.Lms / 2, (self.Lls + self.Nsr * self.Lms), -self.Nsr * self.Lms / 2],
                [-self.Nsr * self.Lms / 2, -self.Nsr * self.Lms / 2, (self.Lls + self.Nsr * self.Lms)]
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

    def IM_getTmec(self):
        '''
        Retorna o torque mecanico do motor de indução
        considerando as 3 fases do motor e os valores
        escalres do Torque de Partida, Nominal e Maximo.

        Retorno:
        Tmec: Vetor do torque mecanico do motor de indução
        Tpartida: Torque de partida do motor de indução (s=1)
        Tnominal: Torque nominal do motor de indução (s=self.s)
        [Tmax, smax]: Torque máximo do motor de indução e o escorregamento correspondente
        '''


        Vs = np.max(self.Vabc[0])
        Tmec = (3 / self.ws) * ((Vs**2) / ((self.Rs + self.Rr / self.sVect) 
                                           + 1j*(self.Xls+self.Xlr))**2) * (self.Rr/self.sVect)
        
        sPartida = 1
        Tpartida = np.abs((3 / self.ws) * ((Vs**2) / ((self.Rs + self.Rr / sPartida) +
                                         1j*(self.Xls+self.Xlr))**2) * (self.Rr / sPartida))

        smax = self.Rr / np.sqrt(self.Rs**2 + (self.Xls+self.Xlr)**2)
        Tmax = np.abs((3 / self.ws) * ((Vs**2) / ((self.Rs + self.Rr / smax) +
                                     1j*(self.Xls+self.Xlr))**2) * (self.Rr / smax))
        

        Tnominal = (3 / self.ws) * ((Vs**2) / ((self.Rs + self.Rr / self.s) 
                                                         + 1j*(self.Xls+self.Xlr))**2) * (self.Rr/self.s)

        return Tmec, Tpartida, np.abs(Tnominal), [Tmax, smax]
    
    def IM_fsSpeedControl(self, f=None):
        '''
        Retorna os torques a partir de uma frequencia

        Parametros:
        f: Vetor de frequencias de controle do motor de indução (Hz)

        Retorno:
        Tmec: Vetor do torque mecanico do motor de indução para cada frequencia de controle
        '''

        if f is None:
            print("Frequencia de controle nao especificada")
            return None

        Vs = np.max(self.Vabc[0])
        ws = 2 * pi * f * (2 / self.Polos)

        Tmec = (3 / ws) * ((Vs**2) / ((self.Rs + self.Rr / self.sVect) 
                                          + 1j*(self.Xls+self.Xlr))**2) * (self.Rr/self.sVect)

        return Tmec

    def IM_getIl(self, HP, Vt, LetraCodigoNominal=None):
        '''
        Retorna a Corrente de partida do motor de inducao a partir
        da Tabela de letras de Código NEMA. 
        
        Baseado no livro "Fundamentos de Maquinas Eletricas; Sthephen J. Chapman
        5 edicao, 2013"
        
        Parametros:
        HP: Potencia do motor em HP
        Vt: Tensão de partida do motor
        LetraCodigoNominal: Letra do codigo nominal do motor (ex: A, B, C, etc)

        Retorno:
        Corrente de partida do motor em Amperes
        '''
        fator = 0

        if LetraCodigoNominal is None:
            print("Codigo Nominal nao especificado")
            return None
        
        match LetraCodigoNominal:
            case 'A':
                fator = 3.15
            case 'B':
                fator = 3.55
            case 'C':
                fator = 4.0
            case 'D':
                fator = 4.5
            case 'E':
                fator = 5.0
            case 'F':
                fator = 5.6
            case 'G':
                fator = 6.3
            case 'H':
                fator = 7.1
            case 'J':
                fator = 8.0
            case 'K':
                fator = 9.0
            case 'L':
                fator = 10.0
            case 'M':
                fator = 11.2
            case 'N':
                fator = 12.5
            case 'P':
                fator = 14.0
            case 'R':
                fator = 16.0
            case 'S':
                fator = 18.0
            case 'T':
                fator = 20.0
            case 'U':
                fator = 22.4
            case 'V':
                fator = 25.0
            case _:
                print("Codigo Nominal nao encontrado")
                return None

        Spartida = HP*fator # kVa

        return (Spartida/(np.sqrt(3)*Vt)*1000)