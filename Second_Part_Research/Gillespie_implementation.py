#%%
from tqdm import tqdm
from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt
#%%
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
    return damage_propensity, Repair_propensity, Mortality_propensity, 

global N_total, tiempo_maximo, a, b, r, s, C, mu

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

        if np.random.rand() > np.abs(maximum_rate - S_T)/maximum_rate:
            
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
N_total = 200
a = 0.03*N_total
b = 0.09
r = 0.8*N_total
s = 1/tiempo_maximo
C = 1.87
mu = 0.0
initial_condition = 0.1

x0 = np.array([0., N_total*initial_condition, 0.])
num_cel = 1000 #número de células 
celulas = np.array([Estado_celula(x0,np.arange(0.,tiempo_maximo,1.)) for i in tqdm(range(num_cel))])

suma = np.nansum(celulas[:, :, 1], axis=0)
longitud_valida = np.sum(~np.isnan(celulas[:, :, 1]), axis=0)
promedio_curva_frailty_index = np.divide(suma, longitud_valida, out=np.zeros_like(suma), where=longitud_valida != 0)

"""
plt.figure(figsize=(8,5))
plt.title(r"Frailty Index Stochastic Simulation", fontsize = 16)
plt.xlabel(r"Time [?]", fontsize = 14)
plt.ylabel(r"Frailty Index", fontsize = 14)
plt.plot(np.arange(0.,tiempo_maximo,1.), promedio_curva_frailty_index/N_total, label = "Average")
plt.legend()
"""

def frailty_index_differential_equation(f, t, a, b, r, s):
    dfdt = a * (1 - f) * (1 + b * t) - f * r * (1 - s * t)
    return dfdt

t = np.linspace(0, tiempo_maximo, 200)  # 100 time steps from 0 to 10
f_solution = odeint(frailty_index_differential_equation, initial_condition, t, args=(a, b, r, s)) #Solution of the differential equation using Odeint

"""
plt.plot(t, f_solution, label='Frailty Index')
plt.xlabel(r'Time')
plt.ylabel(r'Frailty Index')
plt.title(r'Frailty Index')
plt.legend()
"""

plt.figure()
plt.title(r"Frailty Index of population ", fontsize = 16)
plt.plot(promedio_curva_frailty_index/N_total, label = "Gillespie Simulation", color = "red")
plt.plot(t, f_solution, label='ODEs solution', color = "Green")
plt.xlabel(r"Time [Years]", fontsize = 14)
plt.ylabel(r"Frailty Index", fontsize = 14)
plt.legend()
# %%
#Visualización de Frailty Index distribuciones