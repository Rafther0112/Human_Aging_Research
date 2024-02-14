#%%
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import Normalize
from tqdm import tqdm
#%%
moratility_repair_cross_experimentation = np.load('/Users/rafther0112/Documents/GitHub/AGING_RESULTS_SIMULATIONS/Simulacion_cruzado_Damage_Mortality_final_tres_recover.npy', allow_pickle=True)

#%% FRAILTY INDEX DISTRIBUTIONS
frailty_index_curves_General = []
for mortalidad_fija in tqdm(moratility_repair_cross_experimentation):
    frailty_curves_general = []
    for frailty_curves in mortalidad_fija:
        curva_temporal = []
        for tiempo in np.arange(0,80):
            frailty_tiempo_fijo = []
            for persona in np.arange(0,500):
                if frailty_curves[persona, tiempo, 2] == 1:
                    frailty_tiempo_fijo.append(frailty_curves[persona, tiempo, 1])
                else:
                    None
            curva_temporal.append(frailty_tiempo_fijo)
        frailty_curves_general.append(curva_temporal)
    frailty_index_curves_General.append(frailty_curves_general)
#%%
age_of_frailty_index_maximum = []
for mortalidad_fija in tqdm(moratility_repair_cross_experimentation):

    age_curves_general = []
    for frailty_curves in mortalidad_fija:
        age_curva_temporal = []
        for tiempo in np.arange(0,80):
            age_frailty_tiempo_fijo = []
            for persona in np.arange(0,500):
                if frailty_curves[persona, tiempo, 2] == 1:
                    if frailty_curves[persona, tiempo, 0] <= 100 and frailty_curves[persona, tiempo, 0]:
                        age_frailty_tiempo_fijo.append(frailty_curves[persona, tiempo, 0])
                    else:
                        None
                else:
                    None
            age_curva_temporal.append(age_frailty_tiempo_fijo)
        age_curves_general.append(age_curva_temporal)
    age_of_frailty_index_maximum.append(age_curves_general)
# %%

valor_muerte = 0
valor_recovery = 25
for valor_muerte in range(len(valores_de_mu)):

    fig, ax = plt.subplots(figsize = (8,5))

    ax.hist(np.concatenate(frailty_index_curves_General[valor_muerte][valor_recovery]), bins = 15)
    ax.set_title(r"Frailty Index distribution across lifetime" + "\n" + fr"# Nodes: {100} | Recovery: {0.8} | b: {0.09} | s: {0.01} | c: {2.87}" + "\n" + f" Mortality: {round(valores_de_mu[valor_muerte],3)} Recovery: {round(valores_de_recovery[valor_recovery],3)}")
    ax.set_xlabel(r'Frailty Index', fontsize = 14)
    ax.set_ylabel(r'Frecuency', fontsize = 14)
    ax.set_xlim(0,105)
    ax.set_ylim(0,16000)

    plt.savefig(f"Gif_probability_distribution_recovery_mortality/Gif_{valor_muerte}.jpg", dpi = 200)

# %%
valores_de_mu = np.arange(0.5, 1, 0.01)
valores_de_damage = np.arange(0.04, 0.1, 0.001)
percentiles = []
valores_de_mu = np.arange(0.5, 1, 0.01)
for valor_muerte in tqdm(range(len(valores_de_mu))):
    percentile_muerte = []
    for valor_recovery in range(len(valores_de_damage)):
        percentile_muerte.append(np.percentile(np.concatenate(frailty_index_curves_General[valor_muerte][valor_recovery]), 99))
    percentiles.append(percentile_muerte)
percentiles = np.array(percentiles)
#%%
valores_de_mu = np.arange(0.5, 1, 0.01)
valores_de_damage = np.arange(0.04, 0.1, 0.001)
percentiles_Age = []
valores_de_mu = np.arange(0.5, 1, 0.01)
for valor_muerte in tqdm(range(len(valores_de_mu))):
    percentile_muerte = []
    for valor_recovery in range(len(valores_de_damage)):
        percentile_muerte.append(np.percentile(np.concatenate(age_of_frailty_index_maximum[valor_muerte][valor_recovery]),99))
    percentiles_Age.append(percentile_muerte)
percentiles_Age = np.array(percentiles_Age)

#%%
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import Normalize
import numpy as np

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6), sharex=True)

ax1.set_title(r"Frailty Index saturation" + "\n" + fr"# Nodes: {100} | Recovery: {0.3} | b: {0.09} | s: {0.01} | c: {2.87}")
ax1.set_xlabel(r'Damage Rate', fontsize=14)
ax1.set_ylabel(r'Frailty Index Max $(f_{max})$ $[99th]$', fontsize=14)
ax1.set_ylim(0.5, 1)
for i in range(len(percentiles)):
    ax1.plot(valores_de_damage, percentiles[i] / 100, color=plt.cm.Oranges(i / len(percentiles[i])))

norm1 = Normalize(vmin=valores_de_mu[0], vmax=valores_de_mu[-1])  # Adjust the normalization as needed
sm1 = plt.cm.ScalarMappable(cmap=plt.cm.Oranges, norm=norm1)
sm1.set_array([])

divider1 = make_axes_locatable(ax1)
cax1 = divider1.append_axes("right", size="5%", pad=0.1)
cbar1 = plt.colorbar(sm1, cax=cax1)
cbar1.set_label('Mortality Rate')

ax2.set_title(r"Maximum age of death" + "\n" + fr"# Nodes: {100} | Recovery: {0.3} | b: {0.09} | s: {0.01} | c: {2.87}")
ax2.set_xlabel(r'Damage Rate', fontsize=14)
ax2.set_ylabel(r'Time of death max $(\tau_{max})$ $[99th]$', fontsize=14)
ax2.set_ylim(30, 100)
for i in range(len(percentiles)):
    ax2.plot(valores_de_damage, percentiles_Age[i], color=plt.cm.Blues(i / len(percentiles_Age[i])))

norm2 = Normalize(vmin=valores_de_mu[0], vmax=valores_de_mu[-1])  # Adjust the normalization as needed
sm2 = plt.cm.ScalarMappable(cmap=plt.cm.Blues, norm=norm2)
sm2.set_array([])

divider2 = make_axes_locatable(ax2)
cax2 = divider2.append_axes("right", size="5%", pad=0.1)
cbar2 = plt.colorbar(sm2, cax=cax2)
cbar2.set_label('Mortality Rate')
plt.tight_layout()
plt.savefig("Mortality_damage_final.jpg", dpi = 200)
plt.show()

#%%
print(len(percentiles))
# %%
