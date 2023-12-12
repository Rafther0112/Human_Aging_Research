#%%
from tqdm import tqdm
from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt

def damage_Rate(a,b,time,N_individual):
    # return (1- (N_individual/N_total))*a*(1 + b*time)
    return  a*(1 + b*time)

def repair_Rate(r,s,time, N_individual):
    # return (N_individual/N_total)*r*(1 - s*time)
    return r*(1 - s*time)

def mortality_Rate(mu, C, N_individual):
    return mu*(N_individual/N_total)**C

def modelo_constitutivo(a,b,r,s, mu, C, time, N_individual):
    damare_Rate_propensity = damage_Rate(a,b,time, N_individual)
    Repair_Rate_propensity = repair_Rate(r,s,time, N_individual)
    Mortality_Rate_propensity = mortality_Rate(mu, C, N_individual)
    return damare_Rate_propensity, Repair_Rate_propensity, Mortality_Rate_propensity, 

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
        tasa_maxima = a + r 

        τ = (-1/tasa_maxima)*np.log(np.random.rand())
        time+=τ
        if np.random.rand() > np.abs(maximum_rate - S_T)/max_rate:
            
            x = np.random.rand()
            

            if x <= (s_1)/S_T:
                N_individual += 1

            elif x<= (s_1 + s_2)/S_T:
                N_individual -= 1
            
            else: 
                died = True
                return np.array([time, N_individual, died ]) 
        
    return np.array([time, N_individual, died ]) 

def Estado_celula(X0,tiempos):

    X = np.zeros((len(tiempos),len(X0)))
    X[0] = X0
    
    for i in range(1,len(tiempos)):
        X[i] = Gillespie(X[i-1],tiempos[i])
    
    return X
#%%
tiempo_maximo = 1000
N_total = 5000
a = 0.2
b = 0.09
r = 0.8
s = 1/tiempo_maximo
C = 2.87
mu = 1
initial_condition = 0.1
#%%
x0 = np.array([0., N_total*initial_condition, 0.])
num_cel = 1000 #número de células 
celulas = np.array([Estado_celula(x0,np.arange(0.,tiempo_maximo,1.)) for i in tqdm(range(num_cel))])

suma = np.nansum(celulas[:, :, 1], axis=0)
longitud_valida = np.sum(~np.isnan(celulas[:, :, 1]), axis=0)
promedio_curva_frailty_index = np.divide(suma, longitud_valida, out=np.zeros_like(suma), where=longitud_valida != 0)

plt.figure(figsize=(8,5))
plt.title(r"Frailty Index Stochastic Simulation", fontsize = 16)
plt.xlabel(r"Time [?]", fontsize = 14)
plt.ylabel(r"Frailty Index", fontsize = 14)
plt.plot(np.arange(0.,tiempo_maximo,1.), promedio_curva_frailty_index/N_total, label = "Average")
plt.legend()

# %%
s = 0.01
def frailty_index_differential_equation(f, t, a, b, r, s):
    dfdt = a * (1 - f) * (1 + b * t) - f * r * (1 - s * t)
    return dfdt

# Time points
t = np.linspace(0, 100, 200)  # 100 time steps from 0 to 10
f_solution = odeint(frailty_index_differential_equation, initial_condition, t, args=(a, b, r, s)) #Solution of the differential equation using Odeint


plt.figure(figsize=(8,5))
plt.title(r"Frailty Index ODEs solution", fontsize = 16)
plt.xlabel(r"Time [?]", fontsize = 14)
plt.ylabel(r"Frailty Index", fontsize = 14)
plt.plot(t, f_solution, label='Frailty Index')

# %%
