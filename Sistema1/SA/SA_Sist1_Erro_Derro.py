import numpy as np
import matplotlib.pyplot as plt
import skfuzzy as fuzz
from skfuzzy import control
import math
import pickle

#=======================================Função de Discretização do Sistema===========================================================================================================================================================================

def discretizacao(A, B, C, D, ordem, T, n):
    x1 = np.zeros(n)
    x2 = np.zeros(n)
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


    #Sistema de ordem 3:


    for k in range (1,n,1):
        x1[k]=Sig[0,0]*x1[k-1]+ Sig[0,1]*x2[k-1] + B[0,0]*T*u[k-1]
        x2[k]=Sig[1,0]*x1[k-1]+ Sig[1,1]*x2[k-1] + B[1,0]*T*u[k-1]
        y[k]=C[0,0]*x1[k] + C[0,1]*x2[k]
        tempo[k]=tempo[k-1]+T
    return y



#===================================Função Para Calculo do Custo========================================================================================================================================================

def calculaCusto(A, B, C, D, ordem, T, n, erro, derro, saida):
    x1 = np.zeros(n)
    x2 = np.zeros(n)
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
    regra14= control.Rule((~erro['MN']) & (~erro['N']) & (~erro['Z']) & (~erro['P']) & (~erro['MP'])
                          |(~derro['MN']) & (~derro['N']) & (~derro['Z']) & (~derro['P']) & (~derro['MP']), saida['M'])


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
          x1[k]=Sig[0,0]*x1[k-1]+ Sig[0,1]*x2[k-1] + B[0,0]*T*u[k-1]
          x2[k]=Sig[1,0]*x1[k-1]+ Sig[1,1]*x2[k-1] + B[1,0]*T*u[k-1]
          y[k]=C[0,0]*x1[k] + C[0,1]*x2[k]
          tempo[k]=tempo[k-1]+T


#Calculo do IAE

    IAE = np.trapz(np.abs(er), tempo, T)

    return IAE


#=========================================Função de simulação do Controlador Fuzzy======================================================================================================================================================

def simula_controle(erro, derro, saida, y1):
    x1 = np.zeros(n)
    x2 = np.zeros(n)
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
    regra14= control.Rule((~erro['MN']) & (~erro['N']) & (~erro['Z']) & (~erro['P']) & (~erro['MP'])
                          |(~derro['MN']) & (~derro['N']) & (~derro['Z']) & (~derro['P']) & (~derro['MP']), saida['M'])




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
          x1[k]=Sig[0,0]*x1[k-1]+ Sig[0,1]*x2[k-1] + B[0,0]*T*u[k-1]
          x2[k]=Sig[1,0]*x1[k-1]+ Sig[1,1]*x2[k-1] + B[1,0]*T*u[k-1]
          y[k]=C[0,0]*x1[k] + C[0,1]*x2[k]
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
    regra14= control.Rule((~erro['MN']) & (~erro['N']) & (~erro['Z']) & (~erro['P']) & (~erro['MP'])
                          |(~derro['MN']) & (~derro['N']) & (~derro['Z']) & (~derro['P']) & (~derro['MP']), saida['M'])



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
          x1[k]=Sig[0,0]*x1[k-1]+ Sig[0,1]*x2[k-1] + B[0,0]*T*u[k-1]
          x2[k]=Sig[1,0]*x1[k-1]+ Sig[1,1]*x2[k-1] + B[1,0]*T*u[k-1]
          y[k]=C[0,0]*x1[k] + C[0,1]*x2[k]
          tempo[k]=tempo[k-1]+T

    return y, tempo


#==========================Função para criação de funções de pertinência====================================================================================================================================================

def func_pertinencia(pts_erro, pts_derro, pts_saida):

    erro = control.Antecedent(np.arange(-1.0, 1.0, 0.001), 'erro')
    derro = control.Antecedent(np.arange(-0.1, 0.1, 0.0001), 'variacao do erro')
    saida = control.Consequent(np.arange(0,2,0.001), 'saida')

    erro['MN'] = fuzz.trapmf(erro.universe, pts_erro[:4])
    erro['N'] = fuzz.trimf(erro.universe, pts_erro[4:7])
    erro['Z'] = fuzz.trimf(erro.universe, pts_erro[7:10])
    erro['P'] = fuzz.trimf(erro.universe, pts_erro[10:13])
    erro['MP'] = fuzz.trapmf(erro.universe, pts_erro[13:])

    derro['MN'] = fuzz.trimf(derro.universe, pts_derro[:3])
    derro['N'] = fuzz.trimf(derro.universe, pts_derro[3:6])
    derro['Z'] = fuzz.trimf(derro.universe, pts_derro[6:9])
    derro['P'] = fuzz.trimf(derro.universe, pts_derro[9:12])
    derro['MP'] = fuzz.trimf(derro.universe, pts_derro[12:])

    saida['MB'] = fuzz.trapmf(saida.universe, pts_saida[:4])
    saida['B'] = fuzz.trimf(saida.universe, pts_saida[4:7])
    saida['M'] = fuzz.trimf(saida.universe, pts_saida[7:10])
    saida['A'] = fuzz.trimf(saida.universe, pts_saida[10:13])
    saida['MA'] = fuzz.trapmf(saida.universe, pts_saida[13:])

    return erro, derro, saida

#============================================Programa principal(Main)===================================================================================================================================================================



#______________________________________________Modelagem da Planta______________________________________________________________________________________________________________________________________________________________________________

#Planta de controle (n=2)
ordem = 2
A = np.array([(0,1),(-1,-2)])
B = np.array([(0,),(1,)])
C = np.array([(1,0)])
D = 0

T = 0.01  # tempo de amostragem
n = 1000  # numero de amostras

y1 = discretizacao(A, B, C, D, ordem, T, n)

#__________________________________________Parâmetros do Recozimento Simulado________________________________________________________________________________________

maxit = 100                                           # Total de iterações
razao_resfriamento = 0.94                             # Taxa de Resfriamento
temperatura_inicial = 100                             # Temperatura Inicial
delta = 0                                             # Variação do Custo
temperatura_local= temperatura_inicial                # Temperatura Local (Atual)
ner=17                                                # Número de variáveis do erro
nder = 15                                             # Número de variáveis do derro
nsai = 17                                             # Número de variáveis da saida
nvar = 8                                              # Número de variáveis a modificar na solução vizinha
melhor_iter = 0                                       # Melhor iteração


#______________________________________________Inicializando a Solução Inicial________________________________________________________________________________________
   
iga = 0                                         #contador de gerações

pts_erro = np.zeros(ner)           # Vetor com os pontos da função de pertinência erro
pts_derro = np.zeros(nder)         # Vetor com os pontos da função de pertinência variação do erro
pts_saida = np.zeros(nsai)         # Vetor com os pontos da função de pertinência saida

novo_pts_erro = np.zeros(ner)      # Vetor com os pontos da função de pertinência erro vizinho
novo_pts_derro = np.zeros(nder)    # Vetor com os pontos da função de pertinência derro vizinho
novo_pts_saida = np.zeros(nsai)    # Vetor com os pontos da função de pertinência saida vizinho

custo = 0                          # Valor de custo referente ao ISE

#_____________________________________________Inicializando as variaveis Erro, Derro e saida_______________________________________________________________________________________

#Inicializando o Erro

pts_Z = sorted([np.random.uniform(-0.5, 0.3), np.random.uniform(-0.5, 0.5), np.random.uniform(-0.3, 0.5)])
Z_left = pts_Z[0]
Z_center = pts_Z[1]
Z_right = pts_Z[2]

MN_left = -1.0
MN_centerL = -1.0
MN_right = np.random.uniform(-0.9,Z_left)
MN_centerR = np.random.uniform(MN_centerL, MN_right)
pts_MN = (MN_left, MN_centerL, MN_centerR, MN_right)

MP_right = 1.0
MP_centerR = 1.0
MP_left = np.random.uniform(Z_right,0.9)
MP_centerL = np.random.uniform(MP_left, MP_centerR)
pts_MP = (MP_left, MP_centerL, MP_centerR, MP_right)

N_center = np.random.uniform(MN_centerR + 0.01  , Z_center - 0.01)
N_right = max(MN_right,max(N_center, np.random.uniform(Z_left,Z_right)))
N_left = min(Z_left,min(N_center, np.random.uniform(MN_centerR, MN_right)))
pts_N = (N_left, N_center, N_right)

P_center = np.random.uniform(Z_center + 0.01, MP_centerL - 0.01)
P_left = min(MP_left,min(P_center, np.random.uniform(Z_left, Z_right)))
P_right = max(Z_right,max(P_center, np.random.uniform(MP_left,MP_centerL)))
pts_P = (P_left, P_center, P_right)

pts_erro[:4] = pts_MN
pts_erro[4:7] = pts_N
pts_erro[7:10] = pts_Z
pts_erro[10:13] = pts_P
pts_erro[13:] = pts_MP

#Inicializando o Derro

pts_Z = sorted([np.random.uniform(-0.05, 0.03), np.random.uniform(-0.05, 0.05), np.random.uniform(-0.03, 0.05)])
Z_left = pts_Z[0]
Z_center = pts_Z[1]
Z_right = pts_Z[2]

MN_left = -0.1
MN_right = np.random.uniform(-0.09,Z_left)
MN_center = min(np.random.uniform(MN_left, MN_right), Z_center-0.02)
pts_MN = (MN_left, MN_center, MN_right)

MP_right = 0.1
MP_left = np.random.uniform(Z_right,0.09)
MP_center = max(np.random.uniform(MP_left, MP_right), Z_center+0.02)
pts_MP = (MP_left, MP_center, MP_right)

N_center = np.random.uniform(MN_center+0.001 , Z_center-0.001)
N_right = max(MN_right,max(N_center, np.random.uniform(Z_left,Z_center)))
N_left = min(Z_left,min(N_center, np.random.uniform(MN_left, MN_right)))
pts_N = (N_left, N_center, N_right)

P_center = np.random.uniform(Z_center+0.001 , MP_center-0.001)
P_left = min(MP_left,min(P_center, np.random.uniform(Z_center, Z_right)))
P_right = max(Z_right,max(P_center, np.random.uniform(MP_left,MP_right)))
pts_P = (P_left, P_center, P_right)

pts_derro[:3] = pts_MN
pts_derro[3:6] = pts_N
pts_derro[6:9] = pts_Z
pts_derro[9:12] = pts_P
pts_derro[12:] = pts_MP

#Inicializando a Saida

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

pts_saida[:4] = pts_MB
pts_saida[4:7] = pts_B
pts_saida[7:10] = pts_M
pts_saida[10:13] = pts_A
pts_saida[13:] = pts_MA

#__________________________________________________Avaliação da Solução Inicial_____________________________________________________________________________________________

#Criando as funções de pertinência
erro, derro, saida = func_pertinencia(pts_erro, pts_derro, pts_saida)

#Calcula o custo da solução inicial
custo = calculaCusto(A, B, C, D, ordem, T, n, erro, derro, saida)

#Simula o controle fuzzy para a solução inicial
simula_controle(erro, derro, saida, y1)

print('Iteração: 0')
print('Temperatura Local:', temperatura_local)
print('Custo:', custo)

custo_min = np.zeros(maxit)                 # Vetor que contém os valores minimos de custo para cada iteração
custo_min[0] = custo
melhor_custo_min = custo_min[0]             # Armazena o menor custo minimo encontrado
melhor_erro = np.copy(pts_erro)             # Armazena a solução erro que gerou o menor custo
melhor_derro = np.copy(pts_derro)           # Armazena a solução erro que gerou o menor custo
melhor_saida = np.copy(pts_saida)           # Armazena a solução erro que gerou o menor custo
melhor_iter = 0                             #Armazena a iteração que apresento melhor custo

#Atualizando a temperatura
temperatura_local = razao_resfriamento*temperatura_local

#________________________________________________________Começando as Iterações___________________________________________________________________________________________

for iga in range (1,maxit):
#__________________________________________________Calculando uma nova solução vizinha__________________________________________________________________________________________________
    novo_pts_erro = np.copy(melhor_erro)
    novo_pts_derro = np.copy(melhor_derro)
    novo_pts_saida = np.copy(melhor_saida)
    #Solução vizinha erro 
    ncol_erro = np.random.randint(2,15,nvar)                       #escolhe aleatoriamente colunas para serem modificadas
   
   
    for i in range(0,nvar):
        
        if(ncol_erro[i] == 0):
            novo_pts_erro[0]= -1.0
        
        if(ncol_erro[i] == 1):
            novo_pts_erro[1]= -1.0
       
        if(ncol_erro[i] == 2):
            novo_pts_erro[2]= np.random.uniform(novo_pts_erro[1],min(novo_pts_erro[3], novo_pts_erro[5]))
        
        if(ncol_erro[i] == 3):
            novo_pts_erro[3]= np.random.uniform(max(novo_pts_erro[2], novo_pts_erro[4]), novo_pts_erro[6])
        
        if(ncol_erro[i] == 4):
            novo_pts_erro[4]= np.random.uniform(novo_pts_erro[0], min(novo_pts_erro[3], novo_pts_erro[5], novo_pts_erro[7]))
        
        if(ncol_erro[i] == 5):
            novo_pts_erro[5]= np.random.uniform(max(novo_pts_erro[2], novo_pts_erro[4]), min(novo_pts_erro[8], novo_pts_erro[6]))            
        
        if(ncol_erro[i] == 6):
            novo_pts_erro[6]=  np.random.uniform(max(novo_pts_erro[5], novo_pts_erro[7]), novo_pts_erro[9])                
       
        if(ncol_erro[i] == 7):
            novo_pts_erro[7]= np.random.uniform(novo_pts_erro[4], min(novo_pts_erro[6], novo_pts_erro[8], novo_pts_erro[10]))           
        
        if(ncol_erro[i] == 8):
            novo_pts_erro[8]= np.random.uniform(max(novo_pts_erro[5], novo_pts_erro[7]), min(novo_pts_erro[9], novo_pts_erro[11]))
        
        if(ncol_erro[i] == 9):
            novo_pts_erro[9]= np.random.uniform(max(novo_pts_erro[8], novo_pts_erro[10], novo_pts_erro[6]), novo_pts_erro[12])            
        
        if(ncol_erro[i] == 10):
            novo_pts_erro[10]= np.random.uniform(novo_pts_erro[7], min(novo_pts_erro[9], novo_pts_erro[11], novo_pts_erro[13]))
        
        if(ncol_erro[i] == 11):
            novo_pts_erro[11]= np.random.uniform(max(novo_pts_erro[8], novo_pts_erro[10]), min(novo_pts_erro[12], novo_pts_erro[14]))                
        
        if(ncol_erro[i] == 12):
            novo_pts_erro[12]= np.random.uniform(max(novo_pts_erro[11], novo_pts_erro[13], novo_pts_erro[9]), novo_pts_erro[16])                
        
        if(ncol_erro[i] == 13):
            novo_pts_erro[13]= np.random.uniform(novo_pts_erro[10],min(novo_pts_erro[12], novo_pts_erro[14]))                
        
        if(ncol_erro[i] == 14):
            novo_pts_erro[14]= np.random.uniform(max(novo_pts_erro[11], novo_pts_erro[13]), novo_pts_erro[15])                
        
        if(ncol_erro[i] == 15):
            novo_pts_erro[15]= 1.0
                                                                         
        if(ncol_erro[i] == 16):
            novo_pts_erro[16]= 1.0
    
            
    #Solução vizinha erro 
    ncol_derro = np.random.randint(1,14,nvar)                       #escolhe aleatoriamente colunas para serem modificadas

    
    for i in range(0,nvar):
        
        if(ncol_derro[i] == 0):
            novo_pts_derro[0]= -0.1
            
        if(ncol_derro[i] == 1):
            novo_pts_derro[1]= np.random.uniform(novo_pts_derro[0], min(novo_pts_derro[2], novo_pts_derro[4]))
        
        if(ncol_derro[i] == 2):
            novo_pts_derro[2]= np.random.uniform(max(novo_pts_derro[1], novo_pts_derro[3]), novo_pts_derro[5])
        
        if(ncol_derro[i] == 3):
            novo_pts_derro[3]= np.random.uniform(novo_pts_derro[0], min(novo_pts_derro[2], novo_pts_derro[4], novo_pts_derro[6]))
        
        if(ncol_derro[i] == 4):
            novo_pts_derro[4]= np.random.uniform(max(novo_pts_derro[1], novo_pts_derro[3]), min(novo_pts_derro[5], novo_pts_derro[7]))
        
        if(ncol_derro[i] == 5):
            novo_pts_derro[5]= np.random.uniform(max(novo_pts_derro[4],  novo_pts_derro[6], novo_pts_derro[2]), novo_pts_derro[8])
            
        if(ncol_derro[i] == 6):
            novo_pts_derro[6]= np.random.uniform(novo_pts_derro[3], min(novo_pts_derro[5], novo_pts_derro[7], novo_pts_derro[9]))
        
        if(ncol_derro[i] == 7):
            novo_pts_derro[7]= np.random.uniform(max(novo_pts_derro[4], novo_pts_derro[6]), min(novo_pts_derro[8], novo_pts_derro[10]))
        
        if(ncol_derro[i] == 8):
            novo_pts_derro[8]= np.random.uniform(max(novo_pts_derro[7], novo_pts_derro[9], novo_pts_derro[5]), novo_pts_derro[11])
        
        if(ncol_derro[i] == 9):
            novo_pts_derro[9]= np.random.uniform(novo_pts_derro[6], min(novo_pts_derro[8], novo_pts_derro[10], novo_pts_derro[12]))
            
        if(ncol_derro[i] == 10):
            novo_pts_derro[10]= np.random.uniform(max(novo_pts_derro[7], novo_pts_derro[9]), min(novo_pts_derro[11], novo_pts_derro[13]))
            
        if(ncol_derro[i] == 11):
            novo_pts_derro[11]= np.random.uniform(max(novo_pts_derro[10], novo_pts_derro[12], novo_pts_derro[8]), novo_pts_derro[14])
            
        if(ncol_derro[i] == 12):
            novo_pts_derro[12]= np.random.uniform(novo_pts_derro[9], min(novo_pts_derro[11], novo_pts_derro[13]))
            
        if(ncol_derro[i] == 13):
            novo_pts_derro[13]= np.random.uniform(max(novo_pts_derro[10], novo_pts_derro[12]), novo_pts_derro[14])
            
        if(ncol_derro[i] == 14):
            novo_pts_derro[14]= 0.1          
    


#__________________________________________________________Avaliação da nova solução_______________________________________________________________________________________
#Criando as funções de pertinência da solução vizinha
    erro, derro, saida = func_pertinencia(novo_pts_erro, novo_pts_derro, novo_pts_saida)

#Calcula o custo da solução inicial
    novo_custo = calculaCusto(A, B, C, D, ordem, T, n, erro, derro, saida)

#Calculo da variação do custo
    delta = novo_custo - custo
    print('Delta:', delta)

#Aplicando o critério de Metropolis
#Para valores de delta menores que zero a nova solução é melhor que a antiga, logo basta atualizar
    if (delta <= 0):
        pts_erro = novo_pts_erro
        pts_derro = novo_pts_derro
        pts_saida = novo_pts_saida
        custo = novo_custo
        custo_min[iga] = custo
        if(custo < melhor_custo_min):
            melhor_custo_min = custo
            melhor_erro = np.copy(pts_erro)
            melhor_derro = np.copy(pts_derro)
            melhor_saida = np.copy(pts_saida)
            melhor_iter = iga

#Para valores de delta maiores que zero o algoritmo aceita ou não a nova solução atraves de uma probabilidade
    else:
        r = np.random.uniform(0,1)
        p = math.exp(-delta/temperatura_local)
        print('r:',r)
        print('p:', p)

#Probabilidade de aceitação da nova solução
        if (r < p):
            pts_erro = novo_pts_erro
            pts_derro = novo_pts_derro
            pts_saida = novo_pts_saida
            custo = novo_custo
            custo_min[iga] = custo

        else:
            custo_min[iga] = custo

    print('Iteração:',iga)
    print('Temperatura Local:', temperatura_local)
    print('Custo:', custo)
#_______________________________________________Atualização da Temperatura Local_____________________________________________________________________________________

    temperatura_local = razao_resfriamento*temperatura_local

#____________________________________________________Criterio de Parada______________________________________________________________________________________________

    if(melhor_custo_min < 0.1):
      break

#________________________________________________________Resultados__________________________________________________________________________________________________

#Imprime o melhor resultado
print('Melhor valor minimo de custo:', melhor_custo_min)
print('Iteração que apresentou o melhor custo minimo:',melhor_iter)


#Simula o controle fuzzy para a melhor solução encontrada
erro, derro, saida = func_pertinencia(melhor_erro, melhor_derro, melhor_saida)
simula_controle(erro, derro, saida, y1)

#Imprime a evolução do valor do custo durante todas as iterações
plt.plot(custo_min)
plt.xlim(0, maxit)
plt.xlabel('Iterações')
plt.ylabel('Custo (IAE)')
plt.title('SA')
plt.grid(True)
plt.show()


#Simula a curva de controle gerado pelo melhor membro da população
y2, tempo = melhorCurva(erro, derro, saida)

#Calcula o custo do melhor membro da população
melhor_custo = calculaCusto(A, B, C, D, ordem, T, n, erro, derro, saida)
print('Melhor custo:', melhor_custo)


#Impirme a melhor configuração das funções de pertinência encontradas
erro.view()
derro.view()
saida.view()

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

#Comando para armazenar os dados
'''
with open('./novodados2.pkl', 'wb') as arquivo:
    pickle.dump(y2, arquivo)
    
with open('./novocusto2.pkl', 'wb') as arquivo:
    pickle.dump(custo_min, arquivo)
    
'''