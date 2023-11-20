#%%
#Importe de librerias
import numpy as np
from tqdm import tqdm
from numba import jit,njit
import pandas as pd
import json
import matplotlib.pyplot as plt
#___________________________________________________________________________________________________
global total_nodes, tiempo_simulacion
total_nodes = 30
tiempo_simulacion = np.arange(0.,100,1.)

def damage_rate(a, b, t, f):
    R1 = a * (1 - f) * (1 + b * t)
    return R1

def recovery_rate(r, s, t, f):
    R2 = f * r * (1 - s * t)
    return R2

def modelo_constitutivo(time_evolution, frailty_index):

    a = 0.02
    b = 0.09
    r = 0.99
    s = 0.01

    propensidad_damage = damage_rate(a, b, time_evolution, frailty_index)
    propensidad_recovery = recovery_rate(r, s, time_evolution, frailty_index)

    return propensidad_damage, propensidad_recovery


def Gillespie(trp0,tmax):
    """
    Esta funcion se emplea solamente para hacer la evolución de un paso individual en la celula. Evoluciona no un paso temporal, 
    pero si temporalmente la cantidad de veces que pueda evolucionar antes del tmax en una corrida
    """
    t,damage_nodes =trp0 
    frailty = damage_nodes/total_nodes
    while t < tmax:

        s_1, s_2  = modelo_constitutivo(t, frailty)
        S_T = s_1 + s_2 
        
        τ = (-1/S_T)*np.log(np.random.rand())
        x = np.random.rand()

        if x <= (s_1)/S_T:
            damage_nodes += 1
        else:
            damage_nodes -= 1

        t+=τ
    return np.array([t,damage_nodes]) 

def Estado_celula(X0,tiempos):

    X = np.zeros((len(tiempos),len(X0)))
    X[0] = X0
    
    for i in range(1,len(tiempos)):
        X[i] = Gillespie(X[i-1],tiempos[i])
    
    return X

x0 = np.array([20., 2])

num_cel = 1000 #número de células 
celulas = np.array([Estado_celula(x0,tiempo_simulacion) for i in tqdm(range(num_cel))])
celulas_prom = np.mean(celulas,axis=0) #axis = 0 saca el promedio componente a componente de cada célula.

plt.plot(celulas_prom[:,1]/total_nodes)
 

# %%
