No final de cada código de simulação das metaheuristicas existe um campo comentado, que é utilizado para armazenamento dos dados referente as melhor curva de controle encontrada e da curva de minimização.
Para poder armazenar os dados da simulação descomente o codigo e modifique o nome para que o novo arquivo seja salvo com suas especificações (informações do diretorio que quer salva-los. se utilizar somente o nome será salvo na propria pasta do programa).
'''
with open('./dados1.pkl', 'wb') as arquivo:
    pickle.dump(y2, arquivo)
    
with open('./custo1.pkl', 'wb') as arquivo:
    pickle.dump(custo_min, arquivo)
    
'''

No código de comparação dos resultados deve se especificar o caminho e nome do arquivo que será carregado para leitura.

# Carregar dados do arquivo de curvas
with open('./novodados8.pkl', 'rb') as arquivo:
    y1 = pickle.load(arquivo)

Ao especificar somenete o nome ele carregara o arquivo que esta na mesma pasta.