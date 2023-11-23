import numpy as np
import matplotlib.pyplot as plt
import skfuzzy as fuzz
from skfuzzy import control
import math
import pickle
#==========================Função de Discretização do Sistema===========================================================================================================================================================================

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
          

#Calculo do Indice de Desempenho da Integral Absoluta do Erro

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

#Planta de controle (n=2)
ordem = 2
A = np.array([(0,1),(-1,-2)])
B = np.array([(0,),(1,)])
C = np.array([(1,0)])
D = 0

T = 0.01  # tempo de amostragem
n = 1000  # numero de amostras

y1 = discretizacao(A, B, C, D, ordem, T, n)

#____________________________Parâmetros do Algoritmo Genetico___________________________________________________________________________________________________________________________________________________________________________

ner=17                                                # Número de variáveis do erro
nder = 15                                             # Número de variáveis do derro
nsai = 17                                             # Número de variáveis da saida
Ntsai = nsai                                          # Parametro continuo do GA 
Nter =ner                                             # Parametro continuo do GA 
Ntder =nder                                           # Parametro continuo do GA 
popsize= 10                                           # Tamanho da população
mutrate= 0.02                                         # Taxa de mutação
selection= 0.6                                        # Fração da população que será mantida
keep = math.floor(selection*popsize)                  # Membros da população que sobrevivem
nmut = math.ceil((popsize-1)*Nter*mutrate)            # Numero total de mutações 
M = math.ceil((popsize-keep)/2)                       # Numero de cruzamentos

 
#______________________________________Critério de Parada_______________________________________________________________________________________________________________________________________________________________________________

maxit = 100                                     #Número de iterações


#___________________________________Criação das populações iniciais_____________________________________________________________________________________________________________________________________________________________________

iga = 0                                         #contador de gerações
pop_erro = []                                   #população de funções de pertinência do erro
pop_derro = []                                  #população de funções de pertinência da variação do erro
pop_saida = []                                  #população de funções de pertinência da saída

pts_erro = np.zeros((popsize, ner))             #Matriz com os pontos da função de pertinência erro
pts_derro = np.zeros((popsize, nder))           #Matriz com os pontos da função de pertinência variação do erro
pts_saida = np.zeros((popsize, nsai))           #Matriz com os pontos da função de pertinência saida

custo = np.zeros([popsize, 1])                  #Matriz que recebe os valores de custo referente ao ISE
y = np.zeros((maxit,n))

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
   
    
    pts_N = (-0.5, np.random.uniform(-0.4, min(Z_center,0)),  0)
    N_left = pts_N[0]
    M_center = pts_N[1]
    N_right = pts_N[2]
    
    pts_P = (0, np.random.uniform(max(Z_center,0),0.4),  0.5)
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
    
    pts_Z = (-0.05, np.random.uniform(-0.05, 0.05),  0.05)
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
   
    
    pts_N = (-0.1, np.random.uniform(-0.1, min(0,Z_center)),  0)
    N_left = pts_N[0]
    M_center = pts_N[1]
    N_right = pts_N[2]
    
    pts_P = (0, np.random.uniform(max(0,Z_center), 0.1),  0.1)
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
    

  
#________________________________Avaliação da População Inicial_________________________________________________________________________________________________________________________________________________________________________

for i in range (popsize):           #calcula os custos de cada membro da população
   custo[i] = calculaCusto(A, B, C, D, ordem, T, n, pop_erro[i], pop_derro[i], pop_saida[i])    


indices = np.argsort(custo, axis=0)  
custo = np.sort(custo, axis=0)       


#Reorganizando as populações de acordo com seus valores de custo(Ascendente)

a = np.copy(pts_erro)           
b = np.copy(pts_derro)          
c = np.copy(pts_saida)          
j = 0
                   
for i in indices:               
    pts_erro[j] = a[i]
    pts_derro[j] = b[i]
    pts_saida[j] = c[i]
    j+=1

custo_min = np.zeros(maxit)                          #array que contém os valores mínimos de custo para     cada iteração
custo_min[0] = np.min(custo)                         


melhor_custo_min = custo_min[0]                      # guarda o menor custo mínimo encontrado 

melhor_min_pop_erro = np.zeros((popsize, ner))       #Guardar a população erro que gerou o menor custo    
melhor_min_pop_derro = np.zeros((popsize, nder))     #Guardar a população derro que gerou o menor custo    
melhor_min_pop_saida = np.zeros((popsize, nsai))     #Guardar a população saida que gerou o menor custo 


melhor_min_pop_erro = np.copy(pts_erro)                  
melhor_min_pop_derro = np.copy(pts_derro)               
melhor_min_pop_saida = np.copy(pts_saida)             

melhor_geracao_min = 0                                                            
melhor_custo = 0                                        

#Crias as funções de pertinência com os valores dos pontos reorganizados
pop_erro, pop_derro, pop_saida = func_pertinencia(pts_erro, pts_derro, pts_saida)

#Simula o controle fuzzy para população inicial    
for i in range (popsize):
    simula_controle(pop_erro[i], pop_derro[i], pop_saida[i], y1)

print('Geração: 0')
print('Custos:', custo)

#______________________________________Loop de Iterações das Gerações___________________________________________________________________________________________________________________________________________________________________

for iga in range (1,maxit):         
    print('Geração: ', iga)
                                     # Seleção de pares para o cruzamento


    pesos = custo/np.sum(custo)
    pesos = np.reshape(pesos,(-1,1))                   
    pesos= np.flipud(pesos)                                                                           
    chance= np.cumsum(pesos)                       
    chance = np.sort(np.append(chance, 0))                     

    parceiro1=np.random.rand(M)                    
    parceiro2=np.random.rand(M)                    
    ma = np.zeros(M, dtype = 'int')                
    pa = np.zeros(M, dtype = 'int')
    
    #Verifica a aptidão dos membros da população de forma a escolher os que vão realizar o cruzamento
    ic = 0
    while ic < M:
        for id in range(1,popsize):
            if (parceiro1[ic] <= chance[id]) and (parceiro1[ic]>chance[id-1]):
                ma[ic] =id-1
        ic=ic+1
        
    #Verifica a aptidão dos membros da população de forma a escolher os que vão realizar o cruzamento
    ic = 0
    while ic < M:
        for id in range(1,popsize):        
            if (parceiro2[ic]<=chance[id]) and (parceiro2[ic]>chance[id-1]):
                pa[ic]=id-1
        if(pa[ic]==ma[ic]):
            while (pa[ic] == ma[ic]):
                parceiro2[ic]=np.random.rand() 
                for id in range(1,popsize):
                    if (parceiro2[ic]<=chance[id]) and (parceiro2[ic]>chance[id-1]):
                        pa[ic]=id-1
                    
        ic=ic+1
                                        


#_________________________________________Cruzamento aritimético____________________________________________________________________________________________________________________________________________________________________________
   
    #Cruzamento da população Erro
    ic=0
    alpha = np.random.uniform(0, 1)          
    copia_erro = np.copy(pts_erro)           
    for ix in range(popsize-1, keep-1,-2):
        pts_erro[ix,:] = alpha*copia_erro[ma[ic],:] + (1-alpha)*copia_erro[pa[ic],:]      
        if(ix-1>=keep):
            pts_erro[ix-1,:] = alpha*copia_erro[pa[ic],:] + (1-alpha)*copia_erro[ma[ic],:]   
        ic+=1
    
    #Cruzamento da população Derro    
    ic=0
    alpha = np.random.uniform(0, 1)         
    copia_derro = np.copy(pts_derro)        
    for ix in range(popsize-1, keep-1,-2):
        pts_derro[ix,:] = alpha*copia_derro[ma[ic],:] + (1-alpha)*copia_derro[pa[ic],:]      
        if(ix-1>=keep):
            pts_derro[ix-1,:] = alpha*copia_derro[pa[ic],:] + (1-alpha)*copia_derro[ma[ic],:]   
        ic+=1
        
        
#_____________________________________Aplicando mutação na população________________________________________________________________________________________________________________________________________________________________________
    
    #Mutação para população erro
    
    escolhas = [5, 8, 11]
    mut_lin_erro = np.sort(np.random.randint(6,popsize,nmut))           
    mut_col_erro = np.random.choice(escolhas, nmut)                       
   
   
    for i in range(0,nmut):
        
        
        if(mut_col_erro[i] == 5):
            pts_erro[mut_lin_erro[i],mut_col_erro[i]]= np.random.uniform(max(pts_erro[mut_lin_erro[i],2],
                                                                                pts_erro[mut_lin_erro[i],4]),
                                                                            min(pts_erro[mut_lin_erro[i],8],
                                                                                pts_erro[mut_lin_erro[i],6]))             
        if(mut_col_erro[i] == 8):
            pts_erro[mut_lin_erro[i],mut_col_erro[i]]= np.random.uniform(max(pts_erro[mut_lin_erro[i],5],
                                                                                pts_erro[mut_lin_erro[i],7]),
                                                                            min(pts_erro[mut_lin_erro[i],9],
                                                                                pts_erro[mut_lin_erro[i],11]))
        if(mut_col_erro[i] == 11):
            pts_erro[mut_lin_erro[i],mut_col_erro[i]]= np.random.uniform(max(pts_erro[mut_lin_erro[i],8],
                                                                                pts_erro[mut_lin_erro[i],10]),
                                                                            min(pts_erro[mut_lin_erro[i],12],
                                                                                pts_erro[mut_lin_erro[i],14]))                
        
            
    #Mutação para população derro
    escolhas = [4, 7, 10]
    mut_lin_derro = np.sort(np.random.randint(6,popsize,nmut))       
    mut_col_derro = np.random.choice(escolhas, nmut)                  

    
    for i in range(0,nmut):
        
        
        
        if(mut_col_derro[i] == 4):
            pts_derro[mut_lin_derro[i],mut_col_derro[i]]= np.random.uniform(max(pts_derro[mut_lin_derro[i],1],
                                                                                  pts_derro[mut_lin_derro[i],3]),
                                                                              min(pts_derro[mut_lin_derro[i],5],
                                                                                  pts_derro[mut_lin_derro[i],7]))
        
        if(mut_col_derro[i] == 7):
            pts_derro[mut_lin_derro[i],mut_col_derro[i]]= np.random.uniform(max(pts_derro[mut_lin_derro[i],4],
                                                                                  pts_derro[mut_lin_derro[i],6]),
                                                                              min(pts_derro[mut_lin_derro[i],8],
                                                                                  pts_derro[mut_lin_derro[i],10]))
            
        if(mut_col_derro[i] == 10):
            pts_derro[mut_lin_derro[i],mut_col_derro[i]]= np.random.uniform(max(pts_derro[mut_lin_derro[i],7],
                                                                                  pts_derro[mut_lin_derro[i],9]),
                                                                              min(pts_derro[mut_lin_derro[i],11],
                                                                                  pts_derro[mut_lin_derro[i],13]))
            
              
    
    

    #Transforma as populações com os pontos em funções de pertinência para serem usadas no controlador
    pop_erro, pop_derro, pop_saida = func_pertinencia(pts_erro, pts_derro, pts_saida)

#_____________________________Avaliação das População apos o cruzamento e mutações______________________________________________________________________________________________________________________________________________________

    for i in range (popsize):           #calcula os custos de cada membro da população
        custo[i] = calculaCusto(A, B, C, D, ordem, T, n, pop_erro[i], pop_derro[i], pop_saida[i])
    
    indices = np.argsort(custo, axis=0)  
    custo = np.sort(custo, axis=0)       
    
    #Reorganizando as populações de acordo com seus valores de custo(Ascendente)
    
    a = np.copy(pts_erro)           
    b = np.copy(pts_derro)          
    c = np.copy(pts_saida)          
    j = 0
                       
    for i in indices:               
        pts_erro[j] = a[i]
        pts_derro[j] = b[i]
        pts_saida[j] = c[i]
        j+=1
    
    print('Custos:', custo)
#_________________________________Analisando os mínimos para iteração atual_________________________________

    custo_min[iga] = np.min(custo)                

    if(custo_min[iga] < melhor_custo_min):
        melhor_custo_min = custo_min[iga]
        melhor_geracao_min = iga
        melhor_min_pop_erro = np.copy(pts_erro)                  
        melhor_min_pop_derro = np.copy(pts_derro)                
        melhor_min_pop_saida = np.copy(pts_saida)                
    
#___________________________________________Critério de parada____________________________________________  

    if (custo_min[iga]<0.01):
        break

#______________________________________________Resultados Finais____________________________________________

#Imprime os melhores resultados minimos do custo
print('Melhor valor minimo de custo:',melhor_custo_min)
print('Geração que apresentou o melhor custo minimo:',melhor_geracao_min)

#Simula o controle fuzzy para a melhor população encontrada 
pop_erro, pop_derro, pop_saida = func_pertinencia(melhor_min_pop_erro, melhor_min_pop_derro, melhor_min_pop_saida)    

#Imprime a evolução do valor do custo durante todas as iterações
plt.plot(custo_min)
plt.xlim(0, maxit)
plt.xlabel('Iterações')
plt.ylabel('Custo (IAE)')
plt.title('Algoritmo Genético')
plt.grid(True)
plt.show()

#Simula a curva de controle gerado pelo melhor membro da população 
y2, tempo = melhorCurva(pop_erro[0], pop_derro[0], pop_saida[0])

#Calcula o custo do melhor membro da população
melhor_custo = custo[i] = calculaCusto(A, B, C, D, ordem, T, n, pop_erro[0], pop_derro[0], pop_saida[0])
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
#Salva os dados referente a melhor curva encontrada e a evolução doc custo

with open('./dados5, arquivo)

with open('./custo5.pkl', 'wb') as arquivo:
    pickle.dump(custo_min, arquivo)
'''