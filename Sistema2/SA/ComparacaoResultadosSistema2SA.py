import numpy as np
import matplotlib.pyplot as plt
import skfuzzy as fuzz
from skfuzzy import control
import math
import random
import pickle

#==========================Função Para calculo do IAE==========================================================================================================

def calculoIAE(y):
    ref = 1
    er = np.zeros(n)
    
    for k in range (1, n, 1):
          er[k] = ref - y[k-1]
    
    IAE = np.trapz(np.abs(er), tempo, T)
    print('IAE:',np.round(IAE,3))
    
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
    
            
    #Sistema de ordem 3:
        

    for k in range (1,n,1):
        x1[k]=Sig[0,0]*x1[k-1]+ Sig[0,1]*x2[k-1] + Sig[0,2]*x3[k-1] + Sig[0,3]*x4[k-1] + B[0,0]*T*u[k-1] 
        x2[k]=Sig[1,0]*x1[k-1]+ Sig[1,1]*x2[k-1] + Sig[1,2]*x3[k-1] + Sig[1,3]*x4[k-1] + B[1,0]*T*u[k-1] 
        x3[k]=Sig[2,0]*x1[k-1]+ Sig[2,1]*x2[k-1] + Sig[2,2]*x3[k-1] + Sig[2,3]*x4[k-1] + B[2,0]*T*u[k-1]
        x4[k]=Sig[3,0]*x1[k-1]+ Sig[3,1]*x2[k-1] + Sig[3,2]*x3[k-1] + Sig[3,3]*x4[k-1] + B[3,0]*T*u[k-1]
        y[k]=C[0,0]*x1[k] + C[0,1]*x2[k] + C[0,2]*x3[k] + C[0,3]*x4[k]
        tempo[k]=tempo[k-1]+T
    return y

#=======================================Função de calculo do tempo de subida e assentamento======================================================================================================================================================

def calculo_tempos(y,tempo):
    
    y_final = y[-1]  # Último valor da curva de controle
    tempo_subida = 0
    tempo_assentamento = 0
    
    #Calculo do tempo de subida
    y_10 = 0.1*y_final
    y_90 = 0.9*y_final

    # Encontrando o índice do instante de tempo em que a curva de controle atinge 10% do valor final desejado
    ind_10 = np.argmax(y >= y_10)

    # Encontrando o índice do instante de tempo em que a curva de controle atinge 90% do valor final desejado
    ind_90 = np.argmax(y >= y_90)

    # Calculando a diferença de tempo entre esses dois índices, multiplicando pelo tempo de amostragem T
    tempo_subida = (tempo[ind_90] - tempo[ind_10])

    # Calculo do tempo de assentamento
    # Defina a faixa de tolerância
    limite_inferior = y_final-(y_final*0.02)
    limite_superior = y_final+(y_final*0.02)


    # Encontra o índice do instante de tempo em que a curva de controle sai da faixa de tolerância
   
    cont = 0
    for i in range(0,n):
      if(cont == 1):
        break
      if(y[i]>=limite_inferior and y[i]<=limite_superior):
        for j in range (i,n):
          if(not(y[j]>=limite_inferior and y[j]<=limite_superior)):
            break
          if(j==(n-1)):
            tempo_assentamento = tempo[i]
            cont = 1
            break         

    # Calculo da Ultrapassagem Percentual UP%
    y_pico = np.max(y)
    UP = (y_pico - y_final)*100
    UP = np.round(UP,3)
    
    # Imprima os resultados 
    print("Tempo de Subida:", np.round(tempo_subida,3))
    print("Tempo de Assentamento:", np.round(tempo_assentamento,3))
    print("Sobreelevação: ", y_pico)
    print("Ultrapassagem Percentual UP% ={}% ".format(UP))

#____________________________________________________________________________________________________________________________________________________

#Planta de controle (n=4)
ordem = 4
A = np.array([(0, 1, 0, 0), (0, 0, 1, 0), (0, 0, 0, 1), (-64, -120, -70, -15)])
B = np.array([(0,), (0,), (0,), (1,)])
C = np.array([(64, 0, 0, 0)])
D = 0

T = 0.01  # tempo de amostragem
n = 1000  # numero de amostras

maxit= 100

y7 = discretizacao(A, B, C, D, ordem, T, n)
tempo = np.arange(0,10,0.01)

# Carregar dados do arquivo de curvas
with open('./novodados1.pkl', 'rb') as arquivo:
    y1 = pickle.load(arquivo)

with open('./novodados2.pkl', 'rb') as arquivo:
    y2 = pickle.load(arquivo)

with open('./novodados3.pkl', 'rb') as arquivo:
    y3 = pickle.load(arquivo)

with open('./novodados4.pkl', 'rb') as arquivo:
    y4 = pickle.load(arquivo)

with open('./novodados5.pkl', 'rb') as arquivo:
    y5 = pickle.load(arquivo)

with open('./novodados6.pkl', 'rb') as arquivo:
    y6 = pickle.load(arquivo)

# Carregar dados do arquivo de custos
with open('./novocusto1.pkl', 'rb') as arquivo:
    c1 = pickle.load(arquivo)

with open('./novocusto2.pkl', 'rb') as arquivo:
    c2 = pickle.load(arquivo)

with open('./novocusto3.pkl', 'rb') as arquivo:
    c3 = pickle.load(arquivo)

with open('./novocusto4.pkl', 'rb') as arquivo:
    c4 = pickle.load(arquivo)

with open('./novocusto5.pkl', 'rb') as arquivo:
    c5 = pickle.load(arquivo)

with open('./novocusto6.pkl', 'rb') as arquivo:
    c6 = pickle.load(arquivo)

# Definir os limites do eixo x
#limite_inferior = 0
#limite_superior = 6

#Imprime as cusrvas de controle    
plt.figure()

plt.plot(tempo, y1,color='b', label='Cont. Fuzzy 1')
plt.plot(tempo, y2,color='g', label='Cont. Fuzzy 2')
plt.plot(tempo, y3,color='r', label='Cont. Fuzzy 3')
plt.plot(tempo, y4,color='c', label='Cont. Fuzzy 4')
plt.plot(tempo, y5,color='m', label='Cont. Fuzzy 5')
plt.plot(tempo, y6,color='y', label='Cont. Fuzzy 6')
plt.plot(tempo, y7,color='k', label='Malha Aberta')
#plt.xlim(limite_inferior, limite_superior)
plt.legend(loc='best', frameon=False)
plt.xlabel('Tempo (s)')
plt.ylabel('Ganho')
plt.grid(True)
plt.show()

#Imprime os custos por iteração

plt.plot(c1,color='b', label='Cont. Fuzzy 1')
plt.plot(c2,color='g', label='Cont. Fuzzy 2')
plt.plot(c3,color='r', label='Cont. Fuzzy 3')
plt.plot(c4,color='c', label='Cont. Fuzzy 4')
plt.plot(c5,color='m', label='Cont. Fuzzy 5')
plt.plot(c6,color='y', label='Cont. Fuzzy 6')
plt.xlim(0, maxit)
plt.xlabel('Iterações')
plt.ylabel('Custo (IAE)')
plt.legend(loc='best', frameon=False)
plt.grid(True)
plt.show()

#calculo do IAE
calculoIAE(y1)
print("Custo Minimo:",np.round(np.min(c1),3))
calculoIAE(y2)
print("Custo Minimo:",np.round(np.min(c2),3))
calculoIAE(y3)
print("Custo Minimo:",np.round(np.min(c3),3))
calculoIAE(y4)
print("Custo Minimo:",np.round(np.min(c4),3))
calculoIAE(y5)
print("Custo Minimo:",np.round(np.min(c5),3))
calculoIAE(y6)
print("Custo Minimo:",np.round(np.min(c6),3))
calculoIAE(y7)


#Caucula o tempo de assentamento e de subida
calculo_tempos(y1, tempo)
calculo_tempos(y2, tempo)
calculo_tempos(y3, tempo)
calculo_tempos(y4, tempo)
calculo_tempos(y5, tempo)
calculo_tempos(y6, tempo)
calculo_tempos(y7, tempo)