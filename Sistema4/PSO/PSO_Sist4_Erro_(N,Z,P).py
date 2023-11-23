import numpy as np
import matplotlib.pyplot as plt
import skfuzzy as fuzz
from skfuzzy import control
import pickle

#==========================Função de Discretização do Sistema===========================================================================================================================================================================

def discretizacao(A, B, C, D, ordem, T, n):
    x1 = np.zeros(n)
    x2 = np.zeros(n)
    x3 = np.zeros(n)
    x4 = np.zeros(n)
    u =  np.ones(n)
    y = np.zeros(n)
    tempo = np.zeros(n)
    # Metodo para discretizar o sistema pelo Metodo de Euller

    '''
                                    Modelo Discreto:
                                    mat = AT + T
                                    x(k) = matx(k-1) + BTu(k-1)
                                    y(k) = Cx(k) + Du(k)
    '''

    Sig = A * T + np.eye(ordem)



    for k in range (1,n,1):
        x1[k]=Sig[0,0]*x1[k-1]+ Sig[0,1]*x2[k-1] + Sig[0,2]*x3[k-1] + Sig[0,3]*x4[k-1] + B[0,0]*T*u[k-1]
        x2[k]=Sig[1,0]*x1[k-1]+ Sig[1,1]*x2[k-1] + Sig[1,2]*x3[k-1] + Sig[1,3]*x4[k-1] + B[1,0]*T*u[k-1]
        x3[k]=Sig[2,0]*x1[k-1]+ Sig[2,1]*x2[k-1] + Sig[2,2]*x3[k-1] + Sig[2,3]*x4[k-1] + B[2,0]*T*u[k-1]
        x4[k]=Sig[3,0]*x1[k-1]+ Sig[3,1]*x2[k-1] + Sig[3,2]*x3[k-1] + Sig[3,3]*x4[k-1] + B[3,0]*T*u[k-1]
        y[k]=C[0,0]*x1[k] + C[0,1]*x2[k] + C[0,2]*x3[k] + C[0,3]*x4[k]
        tempo[k]=tempo[k-1]+T
    return y



#===================================Função Para Calculo do Custo========================================================================================================================================================

def calculaCusto(A, B, C, D, ordem, T, n, erro, derro, saida):
    x1 = np.zeros(n)
    x2 = np.zeros(n)
    x3 = np.zeros(n)
    x4 = np.zeros(n)
    u =  np.ones(n)
    y = np.zeros(n)
    tempo = np.zeros(n)
    ref = 1
    er = np.zeros(n)


    #Definição das váriavies de entrada e saída do controlador Fuzzy

    erro_crisp = np.zeros(n)
    derro_crisp = np.zeros(n)

    # Metodo para discretizar o sistema pelo Metodo de Euller

    Sig = A * T + np.eye(ordem)



    #Criação de Regras

    regra1  = control.Rule(erro['MN'] & (derro['MP'] | derro['P'])               , saida['MA'])
    regra2  = control.Rule(erro['MN'] & (derro['Z'] | derro['N'] | derro['MN'])  , saida['MB'])
    regra3  = control.Rule(erro['N']  & (derro['MP'] | derro['P'])               , saida['A'])
    regra4  = control.Rule(erro['N']  & (derro['Z'] | derro['N'])                , saida['B'])
    regra5  = control.Rule(erro['N']  &  derro['MN']                             , saida['MB'])
    regra6  = control.Rule(erro['Z']  & derro['MN']                              , saida['B'])
    regra7  = control.Rule(erro['Z']  & (derro['N'] |derro['Z'] | derro['P'])    , saida['M'])
    regra8  = control.Rule(erro['Z']  & derro['MP']                              , saida['A'])
    regra9  = control.Rule(erro['P']  & (derro['Z'] | derro['N'] | derro['MN'])  , saida['A'])
    regra10 = control.Rule(erro['P']  & derro['P']                               , saida['B'])
    regra11 = control.Rule(erro['P']  & derro['MP']                              , saida['MB'])
    regra12 = control.Rule(erro['MP'] & (derro['MP'] | derro['P'])               , saida['MB'])
    regra13 = control.Rule(erro['MP'] & (derro['Z'] | derro['N'] | derro['MN'])  , saida['MA'])
    #Regra de fallback para os casos onde o valor de entrada não pertence a nenhum conjunto fuzzy
    regra14= control.Rule(((~erro['MN']) & (~erro['N']) & (~erro['Z']) & (~erro['P']) & (~erro['MP']))
                          |((~derro['MN']) & (~derro['N']) & (~derro['Z']) & (~derro['P']) & (~derro['MP'])), saida['M'])




    recomendacao_saida = control.ControlSystem([regra1,regra2,regra3,
                                                regra4,regra5,regra6,
                                                regra7,regra8,regra9,
                                                regra10,regra11,regra12,
                                                regra13,regra14])

    controlador = control.ControlSystemSimulation(recomendacao_saida)

    #Laço de simulação do Controlador Fuzzy

    for k in range (1, n, 1):
          er[k] = ref - y[k-1]
          erro_crisp[k] = er[k]
          derro_crisp[k] = 0.1*(erro_crisp[k]-erro_crisp[k-1])
          controlador.input['erro']= erro_crisp[k]
          controlador.input['variacao do erro']= derro_crisp[k]
          controlador.compute()
          u[k] = controlador.output['saida']
          x1[k]=Sig[0,0]*x1[k-1]+ Sig[0,1]*x2[k-1] + Sig[0,2]*x3[k-1] + Sig[0,3]*x4[k-1] + B[0,0]*T*u[k-1]
          x2[k]=Sig[1,0]*x1[k-1]+ Sig[1,1]*x2[k-1] + Sig[1,2]*x3[k-1] + Sig[1,3]*x4[k-1] + B[1,0]*T*u[k-1]
          x3[k]=Sig[2,0]*x1[k-1]+ Sig[2,1]*x2[k-1] + Sig[2,2]*x3[k-1] + Sig[2,3]*x4[k-1] + B[2,0]*T*u[k-1]
          x4[k]=Sig[3,0]*x1[k-1]+ Sig[3,1]*x2[k-1] + Sig[3,2]*x3[k-1] + Sig[3,3]*x4[k-1] + B[3,0]*T*u[k-1]
          y[k]=C[0,0]*x1[k] + C[0,1]*x2[k] + C[0,2]*x3[k] + C[0,3]*x4[k]
          tempo[k]=tempo[k-1]+T


    #Calculo do IAE

    IAE = np.trapz(np.abs(er), tempo, T)

    return IAE


#=========================================Função de simulação do Controlador Fuzzy======================================================================================================================================================

def simula_controle(erro, derro, saida, y1):
    x1 = np.zeros(n)
    x2 = np.zeros(n)
    x3 = np.zeros(n)
    x4 = np.zeros(n)
    u =  np.ones(n)
    y = np.zeros(n)
    tempo = np.zeros(n)
    ref = 1
    er = np.zeros(n)



    #Definição das váriavies de entrada e saída do controlador Fuzzy

    erro_crisp = np.zeros(n)
    derro_crisp = np.zeros(n)

    # Metodo para discretizar o sistema pelo Metodo de Euller

    Sig = A * T + np.eye(ordem)



    #Criação de Regras

    regra1  = control.Rule(erro['MN'] & (derro['MP'] | derro['P'])               , saida['MA'])
    regra2  = control.Rule(erro['MN'] & (derro['Z'] | derro['N'] | derro['MN'])  , saida['MB'])
    regra3  = control.Rule(erro['N']  & (derro['MP'] | derro['P'])               , saida['A'])
    regra4  = control.Rule(erro['N']  & (derro['Z'] | derro['N'])                , saida['B'])
    regra5  = control.Rule(erro['N']  &  derro['MN']                             , saida['MB'])
    regra6  = control.Rule(erro['Z']  & derro['MN']                              , saida['B'])
    regra7  = control.Rule(erro['Z']  & (derro['N'] |derro['Z'] | derro['P'])    , saida['M'])
    regra8  = control.Rule(erro['Z']  & derro['MP']                              , saida['A'])
    regra9  = control.Rule(erro['P']  & (derro['Z'] | derro['N'] | derro['MN'])  , saida['A'])
    regra10 = control.Rule(erro['P']  & derro['P']                               , saida['B'])
    regra11 = control.Rule(erro['P']  & derro['MP']                              , saida['MB'])
    regra12 = control.Rule(erro['MP'] & (derro['MP'] | derro['P'])               , saida['MB'])
    regra13 = control.Rule(erro['MP'] & (derro['Z'] | derro['N'] | derro['MN'])  , saida['MA'])
    #Regra de fallback para os casos onde o valor de entrada não pertence a nenhum conjunto fuzzy
    regra14= control.Rule(((~erro['MN']) & (~erro['N']) & (~erro['Z']) & (~erro['P']) & (~erro['MP']))
                          |((~derro['MN']) & (~derro['N']) & (~derro['Z']) & (~derro['P']) & (~derro['MP'])), saida['M'])




    recomendacao_saida = control.ControlSystem([regra1,regra2,regra3,
                                                regra4,regra5,regra6,
                                                regra7,regra8,regra9,
                                                regra10,regra11,regra12,
                                                regra13,regra14])

    controlador = control.ControlSystemSimulation(recomendacao_saida)

    #Laço de simulação do Controlador Fuzzy


    for k in range (1, n, 1):
          er[k] = ref - y[k-1]
          erro_crisp[k] = er[k]
          derro_crisp[k] = 0.1*(erro_crisp[k]-erro_crisp[k-1])
          controlador.input['erro']= erro_crisp[k]
          controlador.input['variacao do erro']= derro_crisp[k]
          controlador.compute()
          u[k] = controlador.output['saida']
          x1[k]=Sig[0,0]*x1[k-1]+ Sig[0,1]*x2[k-1] + Sig[0,2]*x3[k-1] + Sig[0,3]*x4[k-1] + B[0,0]*T*u[k-1]
          x2[k]=Sig[1,0]*x1[k-1]+ Sig[1,1]*x2[k-1] + Sig[1,2]*x3[k-1] + Sig[1,3]*x4[k-1] + B[1,0]*T*u[k-1]
          x3[k]=Sig[2,0]*x1[k-1]+ Sig[2,1]*x2[k-1] + Sig[2,2]*x3[k-1] + Sig[2,3]*x4[k-1] + B[2,0]*T*u[k-1]
          x4[k]=Sig[3,0]*x1[k-1]+ Sig[3,1]*x2[k-1] + Sig[3,2]*x3[k-1] + Sig[3,3]*x4[k-1] + B[3,0]*T*u[k-1]
          y[k]=C[0,0]*x1[k] + C[0,1]*x2[k] + C[0,2]*x3[k] + C[0,3]*x4[k]
          tempo[k]=tempo[k-1]+T


    erro.view()
    derro.view()
    saida.view()

    plt.figure()
    plt.plot(tempo, y1, label='Malha Aberta')
    plt.plot(tempo, y, label='Controlador Fuzzy')
    plt.title('Resposta do Sistema em Malha Aberta e com Controlador Fuzzy')
    plt.legend(loc='lower right')
    plt.xlabel('Tempo (s)')
    plt.ylabel('Ganho')
    plt.grid(True)
    plt.show()

#=======================================Simula apenas a melhor curva encontrada======================================================================================================================================================

def melhorCurva(erro, derro, saida):
    x1 = np.zeros(n)
    x2 = np.zeros(n)
    x3 = np.zeros(n)
    x4 = np.zeros(n)
    u =  np.ones(n)
    y = np.zeros(n)
    tempo = np.zeros(n)
    ref = 1
    er = np.zeros(n)


    #Definição das váriavies de entrada e saída do controlador Fuzzy

    erro_crisp = np.zeros(n)
    derro_crisp = np.zeros(n)

    # Metodo para discretizar o sistema pelo Metodo de Euller

    Sig = A * T + np.eye(ordem)


    #Criação de Regras

    regra1  = control.Rule(erro['MN'] & (derro['MP'] | derro['P'])               , saida['MA'])
    regra2  = control.Rule(erro['MN'] & (derro['Z'] | derro['N'] | derro['MN'])  , saida['MB'])
    regra3  = control.Rule(erro['N']  & (derro['MP'] | derro['P'])               , saida['A'])
    regra4  = control.Rule(erro['N']  & (derro['Z'] | derro['N'])                , saida['B'])
    regra5  = control.Rule(erro['N']  &  derro['MN']                             , saida['MB'])
    regra6  = control.Rule(erro['Z']  & derro['MN']                              , saida['B'])
    regra7  = control.Rule(erro['Z']  & (derro['N'] |derro['Z'] | derro['P'])    , saida['M'])
    regra8  = control.Rule(erro['Z']  & derro['MP']                              , saida['A'])
    regra9  = control.Rule(erro['P']  & (derro['Z'] | derro['N'] | derro['MN'])  , saida['A'])
    regra10 = control.Rule(erro['P']  & derro['P']                               , saida['B'])
    regra11 = control.Rule(erro['P']  & derro['MP']                              , saida['MB'])
    regra12 = control.Rule(erro['MP'] & (derro['MP'] | derro['P'])               , saida['MB'])
    regra13 = control.Rule(erro['MP'] & (derro['Z'] | derro['N'] | derro['MN'])  , saida['MA'])
    #Regra de fallback para os casos onde o valor de entrada não pertence a nenhum conjunto fuzzy
    regra14= control.Rule(((~erro['MN']) & (~erro['N']) & (~erro['Z']) & (~erro['P']) & (~erro['MP']))
                          |((~derro['MN']) & (~derro['N']) & (~derro['Z']) & (~derro['P']) & (~derro['MP'])), saida['M'])



    recomendacao_saida = control.ControlSystem([regra1,regra2,regra3,
                                                regra4,regra5,regra6,
                                                regra7,regra8,regra9,
                                                regra10,regra11,regra12,
                                                regra13,regra14])

    controlador = control.ControlSystemSimulation(recomendacao_saida)

    #Laço de simulação do Controlador Fuzzy

    for k in range (1, n, 1):
          er[k] = ref - y[k-1]
          erro_crisp[k] = er[k]
          derro_crisp[k] = 0.1*(erro_crisp[k]-erro_crisp[k-1])
          controlador.input['erro']= erro_crisp[k]
          controlador.input['variacao do erro']= derro_crisp[k]
          controlador.compute()
          u[k] = controlador.output['saida']
          x1[k]=Sig[0,0]*x1[k-1]+ Sig[0,1]*x2[k-1] + Sig[0,2]*x3[k-1] + Sig[0,3]*x4[k-1] + B[0,0]*T*u[k-1]
          x2[k]=Sig[1,0]*x1[k-1]+ Sig[1,1]*x2[k-1] + Sig[1,2]*x3[k-1] + Sig[1,3]*x4[k-1] + B[1,0]*T*u[k-1]
          x3[k]=Sig[2,0]*x1[k-1]+ Sig[2,1]*x2[k-1] + Sig[2,2]*x3[k-1] + Sig[2,3]*x4[k-1] + B[2,0]*T*u[k-1]
          x4[k]=Sig[3,0]*x1[k-1]+ Sig[3,1]*x2[k-1] + Sig[3,2]*x3[k-1] + Sig[3,3]*x4[k-1] + B[3,0]*T*u[k-1]
          y[k]=C[0,0]*x1[k] + C[0,1]*x2[k] + C[0,2]*x3[k] + C[0,3]*x4[k]
          tempo[k]=tempo[k-1]+T

    return y, tempo


#==========================Função para criação de funções de pertinência====================================================================================================================================================

def func_pertinencia(pts_erro, pts_derro, pts_saida):

    pop_erro = []
    pop_derro = []
    pop_saida = []


    for i in range(0,popsize):
        erro = control.Antecedent(np.arange(-1.0, 1.0, 0.001), 'erro')
        derro = control.Antecedent(np.arange(-0.1, 0.1, 0.0001), 'variacao do erro')
        saida = control.Consequent(np.arange(0,2,0.001), 'saida')

        erro['MN'] = fuzz.trapmf(erro.universe, pts_erro[i,:4])
        erro['N'] = fuzz.trimf(erro.universe, pts_erro[i,4:7])
        erro['Z'] = fuzz.trimf(erro.universe, pts_erro[i,7:10])
        erro['P'] = fuzz.trimf(erro.universe, pts_erro[i,10:13])
        erro['MP'] = fuzz.trapmf(erro.universe, pts_erro[i,13:])
        pop_erro.append(erro)


        derro['MN'] = fuzz.trimf(derro.universe, pts_derro[i,:3])
        derro['N'] = fuzz.trimf(derro.universe, pts_derro[i,3:6])
        derro['Z'] = fuzz.trimf(derro.universe, pts_derro[i,6:9])
        derro['P'] = fuzz.trimf(derro.universe, pts_derro[i,9:12])
        derro['MP'] = fuzz.trimf(derro.universe, pts_derro[i,12:])
        pop_derro.append(derro)

        saida['MB'] = fuzz.trapmf(saida.universe, pts_saida[i,:4])
        saida['B'] = fuzz.trimf(saida.universe, pts_saida[i,4:7])
        saida['M'] = fuzz.trimf(saida.universe, pts_saida[i,7:10])
        saida['A'] = fuzz.trimf(saida.universe, pts_saida[i,10:13])
        saida['MA'] = fuzz.trapmf(saida.universe, pts_saida[i,13:])
        pop_saida.append(saida)

    return pop_erro, pop_derro, pop_saida

#============================================Programa principal(Main)===================================================================================================================================================================



#______________________________________Modelagem da Planta______________________________________________________________________________________________________________________________________________________________________________

#Planta de controle
ordem = 4

A = np.array([(0,1,0,0),(0,0,1,0),(0,0,0,1),(-5,-106,-121.05,-21.05)])
B = np.array([(0,),(0,),(0,),(1,)])
C = np.array([(55,150,0,0)])
D = 0

T = 0.1  # tempo de amostragem
n = 1000  # numero de amostras

y1 = discretizacao(A, B, C, D, ordem, T, n)

#_____________________________________________________________________________________________________________________
                                    # Inicializando as variáveis

popsize= 10;                                            # tamanho do enxame
ner=17                                                  # Número de variáveis do erro
nder = 15                                               # Número de variáveis do derro
nsai = 17                                               # Número de variáveis da saida
maxit = 200                                             # número máximo de iterações
c1 = 1                                                  # parâmetro cognitivo
c2 = 2                                                  # parâmetro social
melhor_iter = 0                                         # Melhor iteração

#_____________________________________________________________________________________________________________________
                                    # Inicializando as Posições e Velocidades para o Enxame

pop_erro = []                                           #Lista contendo as funções de pertinência da população erro
pop_derro = []                                          #Lista contendo as funções de pertinência da população derro
pop_saida = []                                          #Lista contendo as funções de pertinência da população saida

pts_erro = np.zeros((popsize, ner))                     #Vetor de Posições para a população do Erro
vel_erro = np.zeros((popsize, ner))                     #Vetor de Velocidades para a população do Erro
pts_derro = np.zeros((popsize, nder))                   #Vetor de Posições para a população do Derro
pts_saida = np.zeros((popsize, nsai))                   #Vetor de Posições para a população Saida
pop_global_erro = np.zeros((popsize,ner))               #Matrix que armazena as melhores posições encontradas para o erro

for i in range(popsize):
    erro = control.Antecedent(np.arange(-1.0, 1.0, 0.001), 'erro')

    pts_Z = (-0.25, np.random.uniform(-0.25, 0.25),  0.25)
    Z_left = pts_Z[0]
    Z_center = pts_Z[1]
    Z_right = pts_Z[2]

    pts_MN = (-1,-1,-0.5,-0.25)
    MN_left = pts_MN[0]
    MN_centerL = pts_MN[1]
    MN_right = pts_MN[3]
    MN_centerR = pts_MN[2]


    pts_MP = (0.25,0.5,1,1)
    MP_left = pts_MP[0]
    MP_centerL = pts_MP[1]
    MP_right = pts_MP[3]
    MP_centerR = pts_MP[2]


    pts_N = (-0.5, np.random.uniform(-0.4, min(Z_center-0.1,0)),  0)
    N_left = pts_N[0]
    M_center = pts_N[1]
    N_right = pts_N[2]

    pts_P = (0, np.random.uniform(max(Z_center+0.1,0),0.4),  0.5)
    P_left = pts_P[0]
    P_center = pts_P[1]
    P_right = pts_P[2]

    erro['MN'] = fuzz.trapmf(erro.universe, pts_MN)
    erro['N'] = fuzz.trimf(erro.universe, pts_N)
    erro['Z'] = fuzz.trimf(erro.universe, pts_Z)
    erro['P'] = fuzz.trimf(erro.universe, pts_P)
    erro['MP'] = fuzz.trapmf(erro.universe, pts_MP)

    pop_erro.append(erro)
    pts_erro[i,:4] = pts_MN
    pts_erro[i,4:7] = pts_N
    pts_erro[i,7:10] = pts_Z
    pts_erro[i,10:13] = pts_P
    pts_erro[i,13:] = pts_MP


for i in range(popsize):
    derro = control.Antecedent(np.arange(-0.1, 0.1, 0.0001), 'variacao do erro')

    pts_Z = (-0.05, 0,  0.05)
    Z_left = pts_Z[0]
    Z_center = pts_Z[1]
    Z_right = pts_Z[2]

    pts_MN = (-0.1, -0.1,  -0.05)
    MN_left = pts_MN[0]
    MN_center = pts_MN[1]
    MN_right = pts_MN[2]


    pts_MP = (0.05, 0.1,0.1)
    MP_left = pts_MP[0]
    MP_center = pts_MP[1]
    MP_right = pts_MP[2]


    pts_N = (-0.1, -0.05, 0)
    N_left = pts_N[0]
    M_center = pts_N[1]
    N_right = pts_N[2]

    pts_P = (0, 0.05,  0.1)
    P_left = pts_P[0]
    P_center = pts_P[1]
    P_right = pts_P[2]


    derro['MN'] = fuzz.trimf(derro.universe, pts_MN)
    derro['N'] = fuzz.trimf(derro.universe, pts_N)
    derro['Z'] = fuzz.trimf(derro.universe, pts_Z)
    derro['P'] = fuzz.trimf(derro.universe, pts_P)
    derro['MP'] = fuzz.trimf(derro.universe, pts_MP)

    pop_derro.append(derro)
    pts_derro[i,:3] = pts_MN
    pts_derro[i,3:6] = pts_N
    pts_derro[i,6:9] = pts_Z
    pts_derro[i,9:12] = pts_P
    pts_derro[i,12:] = pts_MP


for i in range(popsize):
    saida = control.Consequent(np.arange(0,2,0.001), 'saida')

    pts_M = (0.75,1,1.25)
    M_left = pts_M[0]
    M_center = pts_M[1]
    M_right = pts_M[2]

    pts_B = (0.5,0.75,1)
    B_center = pts_B[0]
    B_right = pts_B[1]
    B_left = pts_B[2]

    pts_MB = (0,0,0.5,0.75)
    MB_left = pts_MB[0]
    MB_centerL = pts_MB[1]
    MB_right = pts_MB[3]
    MB_centerR = pts_MB[2]

    pts_A = (1,1.25,1.5)
    A_left = pts_A[0]
    A_center = pts_A[1]
    A_right = pts_A[2]


    pts_MA = (1.25,1.5,2,2)
    MA_left = pts_MA[0]
    MA_centerL = pts_MA[1]
    MA_right = pts_MA[3]
    MA_centerR = pts_MA[2]


    saida['MB'] = fuzz.trapmf(saida.universe, pts_MB)
    saida['B'] = fuzz.trimf(saida.universe, pts_B)
    saida['M'] = fuzz.trimf(saida.universe, pts_M)
    saida['A'] = fuzz.trimf(saida.universe, pts_A)
    saida['MA'] = fuzz.trapmf(saida.universe, pts_MA)
    pop_saida.append(saida)

    pts_saida[i, :4] = pts_MB
    pts_saida[i,4:7] = pts_B
    pts_saida[i,7:10] = pts_M
    pts_saida[i,10:13] = pts_A
    pts_saida[i,13:] = pts_MA


#________________________________________________________________________________________________________________________________________________________________________________
                                    # Avaliando a População Inicial

custo = np.zeros([popsize, 1])                #Vetor com os custos de cada membro da população
custo_min = np.zeros(maxit)                   #Vetor que contem os valores minimos de custo
custo_med = np.zeros(maxit)                   #Vetor que contem os valores medios de custo
global_min = np.zeros(maxit)                  #Variavel que armazena o melhor valor de custo ja obtido



#Calculando o custo para as posições iniciais
for i in range (popsize):
    custo[i] = calculaCusto(A, B, C, D, ordem, T, n, pop_erro[i], pop_derro[i], pop_saida[i])


custo_min[0] = np.min(custo)
custo_med[0] = np.mean(custo, axis = 0)
global_min[0] = custo_min[0]

#Crias as funções de pertinência para serem usadas no controlador
pop_erro, pop_derro, pop_saida = func_pertinencia(pts_erro, pts_derro, pts_saida)

#Simula o controle fuzzy para população inicial
for i in range (popsize):
    simula_controle(pop_erro[i], pop_derro[i], pop_saida[i], y1)

print('Iteração: 0')
print('Custo:',custo)


#___________________________________________________________________________________________________________________________________________________________________________________
                                    # Criando as váriaveis locais para atualização das iterações

local_erro = np.copy(pts_erro)                  #localização do minimo local
local_custo = np.copy(custo)                    #Custo da minima local

#____________________________________________________________________________________________________________________________________________________________________________________
                                    # Encontrando a melhor particula na população Inicial

global_custo = np.min(custo)                       #recebe o menor valor de custo encontrado na população inicial
indx = np.where(custo == global_custo)[0][0]       #acha o indice do menor valor de custo no vetor custo

#Inicializa a matriz global de posições com as melhores posições encontradas
for i in range(popsize):
    pop_global_erro[i] = pts_erro[indx]



#______________________________________________________________________________________________________________________________________________________________________________________
                                    # Começando as iterações

for iter in range(1,maxit):
    w = (maxit - iter)/maxit
    r1 = np.random.rand(popsize,ner)
    r2 = np.random.rand(popsize,ner)
    a = np.copy(pts_erro)

#____________________________________Atualização das Velocidades e Posições_____________________________________________________________________________________________________________


    for i in range (popsize):
      for j in range (ner):
        vel_erro[i,j] = (w*vel_erro[i,j] + c1*r1[i,j]*(local_erro[i,j] - pts_erro[i,j]) +c2*r2[i,j]*(pop_global_erro[i,j] - pts_erro[i,j]))

        if (j == 0):
          pts_erro[i,j] = -1.0
          vel_erro[i,j] = pts_erro[i,j] - a[i,j]

        if (j == 1):
          pts_erro[i,j] = -1.0
          vel_erro[i,j] = pts_erro[i,j] - a[i,j]

        if (j == 2):
          pts_erro[i,j] = -0.5
          vel_erro[i,j] = pts_erro[i,j] - a[i,j]

        if (j == 3):
          pts_erro[i,j] = -0.25
          vel_erro[i,j] = pts_erro[i,j] - a[i,j]

        if (j == 4):
          pts_erro[i,j] = -0.5
          vel_erro[i,j] = pts_erro[i,j] - a[i,j]

        if (j == 5):
          pts_erro[i,j] = max(pts_erro[i,2] + 0.1, pts_erro[i, 4], min(pts_erro[i,j]+vel_erro[i,j], pts_erro[i,6]))
          vel_erro[i,j] = pts_erro[i,j] - a[i,j]

        if (j == 6):
          pts_erro[i,j] = 0
          vel_erro[i,j] = pts_erro[i,j] - a[i,j]

        if (j == 7):
          pts_erro[i,j] = -0.25
          vel_erro[i,j] = pts_erro[i,j] - a[i,j]

        if (j == 8):
          pts_erro[i,j] = max(pts_erro[i,5] + 0.1, pts_erro[i, 7], min(pts_erro[i,j]+vel_erro[i,j],pts_erro[i,9]))
          vel_erro[i,j] = pts_erro[i,j] - a[i,j]

        if (j == 9):
          pts_erro[i,j] = 0.25
          vel_erro[i,j] = pts_erro[i,j] - a[i,j]

        if (j == 10):
          pts_erro[i,j] = 0
          vel_erro[i,j] = pts_erro[i,j] - a[i,j]

        if (j == 11):
          pts_erro[i,j] = max(pts_erro[i,8] + 0.1, pts_erro[i, 10], min(pts_erro[i,j]+vel_erro[i,j],pts_erro[i,12]))
          vel_erro[i,j] = pts_erro[i,j] - a[i,j]

        if (j == 12):
          pts_erro[i,j] = 0.5
          vel_erro[i,j] = pts_erro[i,j] - a[i,j]

        if (j == 13):
          pts_erro[i,j] = 0.25
          vel_erro[i,j] = pts_erro[i,j] - a[i,j]

        if (j == 14):
          pts_erro[i,j] = 0.5
          vel_erro[i,j] = pts_erro[i,j] - a[i,j]

        if (j == 15):
          pts_erro[i,j] = 1.0
          vel_erro[i,j] = pts_erro[i,j] - a[i,j]

        if (j == 16):
          pts_erro[i,j] = 1.0
          vel_erro[i,j] = pts_erro[i,j] - a[i,j]


    #Cria as funções de pertinência para serem usadas pelo controlador
    pop_erro, pop_derro, pop_saida = func_pertinencia(pts_erro, pts_derro, pts_saida)


    #Avaliando as novas posições
    for i in range (popsize):           #calcula os custos de cada membro da população e adiciona ao vetor de custo
        custo[i] = calculaCusto(A, B, C, D, ordem, T, n, pop_erro[i], pop_derro[i], pop_saida[i])


#__________________________________________Atualizando as variavés locais_________________________________________________________________________________________________________

    melhor_custo = (custo < local_custo)*1     #Compara os cutos obtidos na iteração anterior com os custos atuais
                                               #e indica (True) ondes os custos atuais são melhores
    pior_custo = (~(custo < local_custo))*1    #Compara os cutos obtidos na iteração anterior com os custos atuais
                                               #e indica (True) ondes os custos atuais são piores



    #Atualizando o custo local
    local_custo = local_custo*pior_custo + custo*melhor_custo      #Atualiza os custos locais com os valores de custo atuais menores


    #Substitui as particulas locais que obtiveram valores de custo piores

    for k in range (popsize):
        if (melhor_custo[k] == 1):
            local_erro[k]=pts_erro[k]



#_______________________________________Atualizando as variaveis globais___________________________________________________________________________________________________________

    custo_min_local = np.min(local_custo)                    #Melhor custoobtido até o momento
    ind = np.where(custo_min_local == local_custo)[0][0]     #Posição no vetor em que o menor custo se encontra

    #Atualização do custo global e posição global
    #Compara o menor custo obtido na iteração atual com o menor custo encontrado anteriormente.
    if (custo_min_local < global_custo):
        for i in range(popsize):
          pop_global_erro[i] = pts_erro[ind]
        indx = ind
        global_custo = custo_min_local
        melhor_iter = iter

    #Menor valor de custo encontrado na Iteração
    custo_min[iter] = np.min(custo)

    #Melhor valor de custo até o momento
    global_min[iter]= global_custo

    #Media do Custo da Iteração
    custo_med[iter] = np.mean(custo)

    #Imprimindo os melhores valores para cada Iteração

    print('iteração:',iter)
    print('Custo:',custo)

#__________________________________________________Critério de parada______________________________________________________________________________________________

    if(custo_min_local < 0.1):
      break


#____________________________________________________Resultados Finais_______________________________________________________________________________________________________________

#Imprime os melhores resultados medios e minimos do custo
print('Melhor valor minimo de custo:', global_custo)
print('Iteração que apresentou o melhor custo minimo:',melhor_iter)


#Imprime a evolução do valor do custo durante todas as iterações
plt.plot(custo_min)
plt.xlim(0, maxit)
plt.xlabel('Iterações')
plt.ylabel('Custo (IAE)')
plt.title('PSO')
plt.grid(True)
plt.show()

#Simula o controle fuzzy para a melhor posição encontrada
pop_erro, pop_derro, pop_saida = func_pertinencia(pop_global_erro, pts_derro, pts_saida)

#Simula a curva de controle gerado pelo melhor membro da população
y2, tempo = melhorCurva(pop_erro[0], pop_derro[0], pop_saida[0])

#Calcula o custo do melhor membro da população
melhor_custo = calculaCusto(A, B, C, D, ordem, T, n, pop_erro[0], pop_derro[0], pop_saida[0])
print('Melhor custo:', melhor_custo)


#Impirme a melhor configuração das funções de pertinência encontradas
pop_erro[0].view()
pop_derro[0].view()
pop_saida[0].view()

#Imprime a meelhor curva de controle encontrada
plt.figure()
plt.plot(tempo, y1, label='Malha Aberta')
plt.plot(tempo, y2, label='Controlador Fuzzy')
plt.title('Resposta do Sistema em Malha Aberta e com Controlador Fuzzy')
plt.legend(loc='lower right')
plt.xlabel('Tempo (s)')
plt.ylabel('Ganho')
plt.grid(True)
plt.show()

'''
with open('./novodados4.pkl', 'wb') as arquivo:
    pickle.dump(y2, arquivo)
    
with open('./novocusto4.pkl', 'wb') as arquivo:
    pickle.dump(custo_min, arquivo)
'''
    
