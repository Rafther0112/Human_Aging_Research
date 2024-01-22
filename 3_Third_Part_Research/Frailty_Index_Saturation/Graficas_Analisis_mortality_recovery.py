#%%
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import Normalize
from tqdm import tqdm
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
from numba import jit,njit,float64,int32
import numba as nb
import pandas as pd
from derivative import dxdt
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable

from derivative import dxdt
from scipy.interpolate import CubicSpline
import numpy as np
import matplotlib.pyplot as plt
#%%
moratility_repair_cross_experimentation = np.load('/Users/rafther0112/Documents/GitHub/AGING_RESULTS_SIMULATIONS/Simulacion_cruzado_Recovery_Mortalidad.npy', allow_pickle=True)

#%%
frailty_index_curves_General = []
for mortalidad_fija in tqdm(moratility_repair_cross_experimentation):
    frailty_curves_general = []
    for frailty_curves in mortalidad_fija:
        curva_temporal = []
        for tiempo in np.arange(0,80):
            valor_frailty_tiempo_especifico = 0
            contador_frailty_especifico = 0
            for persona in np.arange(0,800):
                if frailty_curves[persona, tiempo, 2] != 1:
                    valor_frailty_tiempo_especifico += frailty_curves[persona, tiempo, 1]
                    contador_frailty_especifico += 1
            if  contador_frailty_especifico == 0:
                curva_temporal.append(-1)
            else:
                curva_temporal.append(valor_frailty_tiempo_especifico/contador_frailty_especifico)
        frailty_curves_general.append(curva_temporal)
    frailty_index_curves_General.append(frailty_curves_general)
frailty_index_curves_General = np.array(frailty_index_curves_General)
#%%
valores_de_mu = np.arange(0.000, 1, 0.01)
valores_de_recovery = np.arange(0.000, 1, 0.02)

valores_parametro = np.linspace(0.000, 1, len(valores_de_recovery))
cmap = cm.get_cmap('plasma')

for posicion_mu in range(len(valores_de_mu)):
    posicion_recovery = 0

    fig, ax = plt.subplots(figsize = (8,5))
    ax.set_title(r"Frailty Index with different mortality and recovery rates" + "\n" + fr"# Nodes: {100} | a: {0.05} | b: {0.09} | s: {0.01} | c: {2.87}" + "\n" + f" Mortality: {round(valores_de_mu[posicion_mu],3)}")
    ax.set_xlabel(r'Time', fontsize = 14)
    ax.set_ylabel(r'Frailty Index', fontsize = 14)
    ax.set_ylim(0,1.1)
    ax.axhline(y = 0.7, color = "red")
    #ax.axhline(y = np.max(frailty_index_curves_General[posicion_mu][posicion_recovery]), color = "green")
    ax.set_xlim(18,105)

    for i in range(len(valores_de_recovery)-posicion_recovery):
        color = cmap(valores_parametro[i])

        if -1 in frailty_index_curves_General[posicion_mu][posicion_recovery+i]:
            first_zero_index = np.argmax(frailty_index_curves_General[posicion_mu][posicion_recovery+i] == -1)
            ax.plot(np.arange(20,100,1)[:first_zero_index],frailty_index_curves_General[posicion_mu][posicion_recovery+i][:first_zero_index]/100, color = color, label = f" Mortality: {valores_de_mu[posicion_mu]}" + "\n" + f"Recovery: {round(valores_de_recovery[posicion_recovery+i],3)}")
        else:
            ax.plot(np.arange(20,100,1),frailty_index_curves_General[posicion_mu][posicion_recovery+i]/100, color = color, label = f" Mortality: {valores_de_mu[posicion_mu]}" + "\n" + f"Recovery: {round(valores_de_recovery[posicion_recovery+i],3)}")
    norm = Normalize(vmin=0.000, vmax=1)  # Ajusta los límites del color según tus necesidades
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])  # Necesario para que funcione el colorbar
    cbar = plt.colorbar(sm, ax=ax, label=r'Recovery Rate $r$')

    plt.savefig(f"Gif_results_Mortality_recovery/Gif_{posicion_mu}.jpg", dpi = 200)

#%%
mortality_rate_general = []
for mortality in tqdm(moratility_repair_cross_experimentation):
    mortalidad_mortality = []
    for repair in mortality:
        mortalidad_repair = []
        for i in np.arange(0,80):
            mortalidad_repair.append(np.sum(repair[:,i,2]))
        mortalidad_mortality.append(mortalidad_repair)
    mortality_rate_general.append(mortalidad_mortality)
mortality_rate_general = np.array(mortality_rate_general)
#%%
plt.plot(mortality_rate_general[5][20]/800)

#%%
valores_de_mu = np.arange(0.000, 1, 0.01)
valores_de_recovery = np.arange(0.000, 1, 0.02)

valores_parametro = np.linspace(0.000, 1, len(valores_de_recovery))
cmap = cm.get_cmap('plasma')

for posicion_mu in range(len(valores_de_mu)):
    posicion_recovery = 0

    fig, ax = plt.subplots(figsize = (8,5))
    ax.set_title(r"Mortality with different mortality and recovery rates" + "\n" + fr"# Nodes: {100} | a: {0.05} | b: {0.09} | s: {0.01} | c: {2.87}" + "\n" + f" Mortality: {round(valores_de_mu[posicion_mu],3)}")
    ax.set_xlabel(r'Time', fontsize = 14)
    ax.set_ylabel(r'Mortality', fontsize = 14)
    ax.set_ylim(0,1.1)
    ax.axhline(y = 0.7, color = "red")
    #ax.axhline(y = np.max(frailty_index_curves_General[posicion_mu][posicion_recovery]), color = "green")
    ax.set_xlim(18,105)

    for i in range(len(valores_de_recovery)-posicion_recovery):
        color = cmap(valores_parametro[i])

        ax.plot(np.arange(20,100,1),mortality_rate_general[posicion_mu][i]/800, color = color)

    norm = Normalize(vmin=0.000, vmax=1)  # Ajusta los límites del color según tus necesidades
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])  # Necesario para que funcione el colorbar
    cbar = plt.colorbar(sm, ax=ax, label=r'Recovery Rate $r$')

    plt.savefig(f"Gif_results_mortality_with_diff_mortality_rate/Gif_{posicion_mu}.jpg", dpi = 200)
# %%
len(mortality_rate_general[10][0][20:40])


#%%
beta_values

#%%
mortality_rate_general.shape
#%%
valores_de_mu = np.arange(0.000, 1, 0.01)

R_parameter_gompertz_law_40_60 = np.zeros((len(valores_de_mu),len(valores_de_recovery)), dtype=object)
R_parameter_gompertz_law_60_80 = np.zeros((len(valores_de_mu),len(valores_de_recovery)), dtype=object)

betas_gompertz_law_40_60 = np.zeros((len(valores_de_mu),len(valores_de_recovery)), dtype=object)
betas_gompertz_law_60_80 = np.zeros((len(valores_de_mu),len(valores_de_recovery)), dtype=object)

for posicion_mortality in tqdm(range(len(valores_de_mu))):
    for posicion_recovery in range(len(valores_de_recovery)):
        betas_40_60 = []
        R_values_40_60 = []
        edad_40_60 = np.arange(40,60,1)
        data_40_60 = mortality_rate_general[posicion_mortality][posicion_recovery][20:40]
        Cubic_data = CubicSpline(edad_40_60, data_40_60, bc_type='natural')
        Spline_data_40_60 = Cubic_data(edad_40_60)
        

        derivative_process = dxdt(Spline_data_40_60, edad_40_60, kind="finite_difference", k=1)
        beta_values = derivative_process/Spline_data_40_60
        R_values = Spline_data_40_60/np.e**(beta_values*edad_40_60)

        betas_gompertz_law_40_60[posicion_mortality][posicion_recovery] = np.array(beta_values)
        R_parameter_gompertz_law_40_60[posicion_mortality][posicion_recovery] = np.array(R_values)


        betas_60_80 = []
        R_values_60_80 = []
        edad_60_80 = np.arange(60,80,1)
        data_60_80 =  mortality_rate_general[posicion_mortality][posicion_recovery][40:60]
        Cubic_data = CubicSpline(edad_60_80, data_60_80, bc_type='natural')
        Spline_data_60_80 = Cubic_data(edad_60_80)
        
        derivative_process = dxdt(Spline_data_60_80, edad_60_80, kind="finite_difference", k=1)
        beta_values = derivative_process/Spline_data_60_80
        R_values = Spline_data_60_80/np.e**(beta_values*edad_60_80)

        betas_gompertz_law_60_80[posicion_mortality][posicion_recovery] = beta_values
        R_parameter_gompertz_law_60_80[posicion_mortality][posicion_recovery] = R_values
#%%
posicion_mortality = 0
posicion_recovery = 35


for i in range(len(valores_de_mu)-posicion_mortality):

    plt.figure(figsize = (8,5))
    plt.title(r"SM Correlation" + "\n" + fr"# Nodes: {100} | a: {0.05} | b: {0.09} | s: {0.01} | c: {2.87}" +"\n"+rf"Mortality: {valores_de_mu[i]} | Recovery: {valores_de_recovery[posicion_recovery]}")
    plt.yscale("log")
    plt.ylim(10**-2.5, 10**3.2)
    plt.xlim(-0.01,0.2)
    plt.xlabel(r"$\alpha$ parameter", fontsize = 14)
    plt.ylabel(r"$Ln(R_0)$ parameter", fontsize = 14)
    plt.scatter(betas_gompertz_law_40_60[i][posicion_recovery], R_parameter_gompertz_law_40_60[i][posicion_recovery], color="red",label = "Ages Window: 40-60")
    plt.scatter(betas_gompertz_law_60_80[i][posicion_recovery], R_parameter_gompertz_law_60_80[i][posicion_recovery], color= "blue", marker='*',label = "Ages Window: 60-80")
    plt.legend()

    plt.savefig(f"SM_correlation/Recovery_{posicion_recovery}/Gif_{i}.jpg", dpi = 200)

#%%
plt.scatter(betas_40_60, R_values_40_60)
plt.scatter(betas_60_80, R_values_60_80)
plt.yscale("log")
#%%

R_parameter_gompertz_law_40_60.append(R_values_40_60)
R_parameter_gompertz_law_60_80.append(R_values_60_80)

betas_gompertz_law_40_60.append(betas_40_60)
betas_gompertz_law_60_80.append(betas_60_80)



#%%
#%%
plt.figure()
plt.title("Frailty Index Distributions")
plt.hist(moratility_repair_cross_experimentation[45][-1][:,40,1], density = True, bins = 10)
plt.xlim(0,100)
# %%

moratility_repair_cross_experimentation[45][-1][:,40,1]
# %%
data = moratility_repair_cross_experimentation[45][-1][:,40,1]/100

