""" 
Here we will try to explain the main reason of the frailty index maximum limit. Therefore we will change the mortality rate in the stochastic process
We will try with the explicit model of damage nodes

tiempo_maximo = 100
N_total = 100
a = 0.02*N_total
b = 0.09
r = 0.9*N_total
s = (1/tiempo_maximo)
C = 3
initial_condition = 0.04

Vamos a modificar mu 
"""
#%%
from tqdm import tqdm
from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt
#%%
global N_total, tiempo_maximo, d, alpha, r, s, C, mu
def damage_Rate(d,alpha,frailty):
    """
    Args:
        a (_type_): _description_
        b (_type_): _description_
        time (_type_): Temporal evolution of the simulation

    Returns:
        _type_: Damage rate at which deterioration of system nodes occurs.
    """
    return  (d + alpha*frailty)
def repair_Rate(r,s,frailty):
    """

    Args:
        r (_type_): _description_
        s (_type_): _description_
        time (_type_): Temporal evolution of the simulation

    Returns:
        _type_: Repair rate at which repair of system nodes occurs.
    """
    return r*(1 - s*frailty)
def mortality_Rate(mu, C, N_individual):
    """

    Args:
        mu (_type_): _description_
        C (_type_): _description_
        N_individual (_type_): _description_

    Returns:
        _type_: _description_
    """
    return mu*(N_individual/N_total)**C
def modelo_constitutivo(d,alpha,r,s, mu, C, N_individual):
    damage_propensity = ((N_total - N_individual))*damage_Rate(d,alpha,N_individual/N_total)
    Repair_propensity = ((N_individual))*repair_Rate(r,s,N_individual/N_total)
    Mortality_propensity = mortality_Rate(mu, C, N_individual)
    return damage_propensity, Repair_propensity, Mortality_propensity
def Gillespie(trp0,tmax):
    """
    Esta funcion se emplea solamente para hacer la evolución de un paso individual en el individuo. Evoluciona no un paso temporal, 
    pero si temporalmente la cantidad de veces que pueda evolucionar antes del tmax en una corrida
    """
    time , N_individual, died  =trp0 

    while time < tmax and not died:
        s_1, s_2, s_3  =  modelo_constitutivo(d,alpha,r,s, mu, C, N_individual)
        S_T = s_1 + s_2 + s_3 

        #maximum_rate = ((N_total - N_individual))*damage_Rate(a,b,tiempo_maximo) + ((N_individual))*repair_Rate(r,s,0) + s_3 
        τ = (-1/S_T)*np.log(np.random.rand())
        time+=τ

        if time < tmax:
            
            x = np.random.rand()
        
            if x <= (s_1)/S_T:
                N_individual += 1

            elif x<= (s_1 + s_2)/S_T:
                N_individual -= 1
            
            else: 
                died = True
                return np.array([time, N_individual, died ]) 
        else: 
            None
    return np.array([time, N_individual, died ]) 
def Estado_celula(X0,tiempos):
    X = np.zeros((len(tiempos),len(X0)))
    X[0] = X0
    
    for i in range(1,len(tiempos)):
        X[i] = Gillespie(X[i-1],tiempos[i])
    
    return X
#%%
tiempo_maximo = 100
N_total = 100
alpha = 0.05
s = (1/tiempo_maximo)
C = 2.87 
mu = 1.48

initial_condition = 0.045
x0 = np.array([20.44, int(N_total*initial_condition), 0.])
num_cel = 1000 #número de células 


valores_de_damage = np.arange(0.00001, 0.001, 0.0001)
valores_de_recovery = np.arange(0.001, 0.01, 0.001)
array_principal = np.zeros((len(valores_de_damage),) + (len(valores_de_recovery),) + (num_cel,80,3 ))

#%%
for posicion_damage, d in enumerate(tqdm(valores_de_damage)):
    for posicion_recovery, r in enumerate(tqdm((valores_de_recovery))):
        array_principal[posicion_damage][posicion_recovery] = np.array([Estado_celula(x0,np.arange(20.,tiempo_maximo,1.)) for i in (range(num_cel))])
#        np.save('/Users/rafther0112/Documents/GitHub/AGING_RESULTS_SIMULATIONS/Simulacion_cruzado_Damage_Mortality_PowerLaw_3.npy', array_principal)
# %%

