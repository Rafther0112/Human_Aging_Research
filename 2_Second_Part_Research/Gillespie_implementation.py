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
    damage_propensity = ((N_total - N_individual))*damage_Rate(a,b,time)
    Repair_propensity = ((N_individual))*repair_Rate(r,s,time)
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
        maximum_rate = ((N_total - N_individual))*damage_Rate(a,b,tiempo_maximo) + ((N_individual))*repair_Rate(r,s,0) + s_3 
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
b = 0.09
s = (1/tiempo_maximo)
C = 2.87
mu = 0.
a = 0.005
initial_condition = 0.045
x0 = np.array([20.44, int(N_total*initial_condition), 0.])
num_cel = 800 #número de células 

valores_de_r = np.linspace(0.1, 0.99, 101)
array_principal = np.zeros((len(valores_de_r),) + (num_cel,tiempo_maximo,3 ))

r = 0.8

celulas = np.array([Estado_celula(x0,np.arange(20.,tiempo_maximo,1.)) for i in tqdm(range(num_cel))])

suma = np.nansum(celulas[:, :, 1], axis=0)
longitud_valida = np.sum(~np.isnan(celulas[:, :, 1]), axis=0)
promedio_curva_frailty_index = np.divide(suma, longitud_valida, out=np.zeros_like(suma), where=longitud_valida != 0)


def frailty_index_differential_equation(f, t, a, b, r, s):
    dfdt = a * (1 - f) * (1 + b * t) - f * r * (1 - s * t)
    return dfdt

t = np.linspace(20, tiempo_maximo, 201)  # 100 time steps from 0 to 10
f_solution = odeint(frailty_index_differential_equation, initial_condition, t, args=(a, b, r, s)) #Solution of the differential equation using Odeint


plt.figure()
plt.title(r"Frailty Index of population ", fontsize = 16)
plt.plot(np.arange(20.,tiempo_maximo,1.), promedio_curva_frailty_index/N_total, label = rf"GS N = {N_total} , $\mu_0$ = 0", color = "red")
plt.plot(t, f_solution, label='ODEs solution')
plt.xlabel(r"Time [Years]", fontsize = 14)
plt.ylabel(r"Frailty Index", fontsize = 14)
plt.legend()
#plt.savefig("Differents_Values_Mortality.jpg", dpi = 1000)

# %%
mortality_data = []
for i in np.arange(0,100):
    mortality_data.append(np.sum(celulas[:,i,2]))
mortality_data = np.array(mortality_data)
# %%

# %%
def frailty_index_differential_equation(f, t, a, b, r, s):
    dfdt = (1 - f)*a*(1 + b*t) - f*r*(1 - s*t)
    return dfdt

t = np.linspace(0, tiempo_maximo, 201)  # 100 time steps from 0 to 10
f_solution = odeint(frailty_index_differential_equation, initial_condition, t, args=(a, b, r, s)) #Solution of the differential equation using Odeint
mu = 0.5
gompertz_law = mu*(f_solution**C)
#%%

plt.figure()
plt.title(r"Mortality Rate  of population ", fontsize = 16)
plt.plot(mortality_data[0:-1]/num_cel, label = rf"GS N = {N_total} , $\mu_0$ = 0.01", color = "crimson")

plt.xlabel(r"Time [Years]", fontsize = 14)
plt.ylabel(r"Mortality Rate", fontsize = 14)
plt.legend()
plt.savefig("Mortality_Curve_scolanada.jpg", dpi = 1000)
# %%

# %%
t = np.linspace(0, tiempo_maximo, 201)
# %%
t
# %%
