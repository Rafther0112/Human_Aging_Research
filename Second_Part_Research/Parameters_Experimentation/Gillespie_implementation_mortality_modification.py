""" 
Aqui vamos a hacer la experimentación para modificar el coeficiente de mortalidad mu y mirar la variación que presenta el frailty index con respecto a la
solución del sistema de ecuaciones diferenciales no acoplado con la mortalidad. 
Vamos a evidenciar que hay una tendencia a los datos aplanarse antes del tiempo cuando se aumenta el coeficiente de mortalidad. 

tiempo_maximo = 100
N_total = 100
a = 0.02*N_total
b = 0.09
r = 0.9*N_total
s = (1/tiempo_maximo)
C = 2.87
initial_condition = 0.04

Vamos a modificar mu 
"""
#%%
from tqdm import tqdm
from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt
#%%
global N_total, tiempo_maximo, a, b, r, s, C, mu
def damage_Rate(a,b,time):
    """
    Args:
        a (_type_): _description_
        b (_type_): _description_
        time (_type_): Temporal evolution of the simulation

    Returns:
        _type_: Damage rate at which deterioration of system nodes occurs.
    """
    return  a*(1 + b*time)
def repair_Rate(r,s,time):
    """

    Args:
        r (_type_): _description_
        s (_type_): _description_
        time (_type_): Temporal evolution of the simulation

    Returns:
        _type_: Repair rate at which repair of system nodes occurs.
    """
    return r*(1 - s*time)
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
def modelo_constitutivo(a,b,r,s, mu, C, time, N_individual):
    damage_propensity = ((N_total - N_individual)/(N_total))*damage_Rate(a,b,time)
    Repair_propensity = ((N_individual)/(N_total))*repair_Rate(r,s,time)
    Mortality_propensity = mortality_Rate(mu, C, N_individual)
    return damage_propensity, Repair_propensity, Mortality_propensity
def Gillespie(trp0,tmax):
    """
    Esta funcion se emplea solamente para hacer la evolución de un paso individual en el individuo. Evoluciona no un paso temporal, 
    pero si temporalmente la cantidad de veces que pueda evolucionar antes del tmax en una corrida
    """
    time , N_individual, died  =trp0 

    while time < tmax and not died:
        s_1, s_2, s_3  =  modelo_constitutivo(a,b,r,s, mu, C, time, N_individual)
        S_T = s_1 + s_2 + s_3 
        maximum_rate = ((N_total - N_individual)/(N_total))*damage_Rate(a,b,tiempo_maximo) + ((N_individual)/(N_total))*repair_Rate(r,s,0) + s_3 
        τ = (-1/maximum_rate)*np.log(np.random.rand())
        time+=τ

        if time < tmax and np.random.rand() > np.abs(maximum_rate - S_T)/maximum_rate:
            
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
a = 0.05*N_total
b = 0.09
r = 0.9*N_total
s = (1/tiempo_maximo)
C = 2.87
initial_condition = 0.04
x0 = np.array([0., N_total*initial_condition, 0.])
num_cel = 1000 #número de células 

valores_de_mu = np.linspace(0.01, 0.9, 90)
array_principal = np.zeros((len(valores_de_mu),) + (num_cel,tiempo_maximo,3 ))

#%%
for posicion_mu, mu in enumerate(tqdm(valores_de_mu)):
    array_principal[posicion_mu] = np.array([Estado_celula(x0,np.arange(0.,tiempo_maximo,1.)) for i in (range(num_cel))])
#%%
np.save('/Users/rafther0112/Documents/GitHub/AGING_RESULTS_SIMULATIONS/Simulacion_modificacion_tasa_mortalidad.npy', array_principal)

