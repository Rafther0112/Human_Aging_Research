#%%
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import Normalize
from tqdm import tqdm


#%%
damage_repair_cross_experimentation = np.load('/Users/rafther0112/Documents/GitHub/AGING_RESULTS_SIMULATIONS/Simulacion_cruzado_Damage_Recovery_conmu_tres.npy', allow_pickle=True)
# %%
frailty_index_curves_General = []
for mortalidad_fija in tqdm(damage_repair_cross_experimentation):
    frailty_curves_general = []
    for frailty_curves in mortalidad_fija:
        curva_temporal = []
        for tiempo in np.arange(0,80):
            valor_frailty_tiempo_especifico = 0
            contador_frailty_especifico = 0
            for persona in np.arange(0,200):
                if frailty_curves[persona, tiempo, 2] != 1:
                    valor_frailty_tiempo_especifico += frailty_curves[persona, tiempo, 1]
                    contador_frailty_especifico += 1
            if  contador_frailty_especifico == 0:
                curva_temporal.append(-1)
            else:
                curva_temporal.append(valor_frailty_tiempo_especifico/contador_frailty_especifico)
        frailty_curves_general.append(curva_temporal)
    frailty_index_curves_General.append(frailty_curves_general)
frailty_index_curves_General = np.array(frailty_index_curves_General)/100
#%%
valores_de_damage = np.arange(0.000, 0.05, 0.001)

valores_de_mu = np.arange(0.000, 1, 0.01)
valores_de_recovery = np.arange(0.000, 1, 0.02)

valores_parametro = np.linspace(0.000, 1, len(valores_de_recovery))
cmap = cm.get_cmap('plasma')

for posicion_damage in range(len(valores_de_damage)):

    posicion_recovery = 0

    fig, ax = plt.subplots(figsize = (8,5))
    ax.set_title(r"Frailty Index with different mortality and recovery rates" + "\n" + fr"# Nodes: {100} | a: {0.05} | b: {0.09} | s: {0.01} | c: {2.87} | mu: {0.03}" + "\n" + f" Damage: {round(valores_de_damage[posicion_damage],3)}")
    ax.set_xlabel(r'Time', fontsize = 14)
    ax.set_ylabel(r'Frailty Index', fontsize = 14)
    ax.set_ylim(0,1.1)
    ax.axhline(y = 0.7, color = "red")
    ax.axhline(y = np.max(frailty_index_curves_General[posicion_damage][posicion_recovery]), color = "green")


    for i in range(len(valores_de_recovery)-posicion_recovery):
        color = cmap(valores_parametro[i])
        ax.plot(np.arange(20,100,1),frailty_index_curves_General[posicion_damage][posicion_recovery+i], color = color, label = f" Mortality: {valores_de_mu[posicion_damage]}" + "\n" + f"Recovery: {round(valores_de_recovery[posicion_recovery+i],3)}")
    norm = Normalize(vmin=0.000, vmax=1)  # Ajusta los límites del color según tus necesidades
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])  # Necesario para que funcione el colorbar
    cbar = plt.colorbar(sm, ax=ax, label=r'Recovery Rate $r$')

    plt.savefig(f"Gif_results_damage_recovery_dos/Gif_{posicion_damage}.jpg", dpi = 200)


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
posicion_mu = 20
posicion_recovery = -1


plt.figure(figsize=(8,5))
plt.title(r"Mortality with different mortality and recovery rates" + "\n" + fr"# Nodes: {100} | a: {0.05} | b: {0.09} | s: {0.01} | c: {2.87}")
plt.ylim(-0.05,1.01)
#plt.axhline(y = 0.7, color = "red")
#plt.axhline(y = np.max(mortality_rate_general[posicion_mu][posicion_recovery]), color = "green")
plt.ylabel(r"Frailty Index FI", fontsize = 16)
plt.xlabel(r"Time [Years]", fontsize = 16)
plt.scatter(np.arange(20,100,1),mortality_rate_general[posicion_mu][posicion_recovery]/800, label = f" Mortality: {valores_de_mu[posicion_mu]}" + "\n" + f"Recovery: {round(valores_de_recovery[posicion_recovery],3)}")
plt.legend()
# %%
plt.figure()
plt.title("Frailty Index Distributions")
plt.hist(moratility_repair_cross_experimentation[45][-1][:,40,1], density = True, bins = 10)
plt.xlim(0,100)
# %%

moratility_repair_cross_experimentation[45][-1][:,40,1]
# %%
data = moratility_repair_cross_experimentation[45][-1][:,40,1]/100

