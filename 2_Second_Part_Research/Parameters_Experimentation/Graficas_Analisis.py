#%%
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import Normalize
from tqdm import tqdm
#%%
mortality_rate_modification = np.load('/Users/rafther0112/Documents/GitHub/AGING_RESULTS_SIMULATIONS/Simulacion_modificacion_tasa_mortalidad.npy', allow_pickle=True)
#%%
frailty_curves_general = []
for frailty_curves in tqdm(mortality_rate_modification):
    curva_temporal = []
    for tiempo in np.arange(0,100):
        valor_frailty_tiempo_especifico = 0
        contador_frailty_especifico = 0
        for persona in np.arange(0,3000):
            if frailty_curves[persona, tiempo, 2] != 1:
                valor_frailty_tiempo_especifico += frailty_curves[persona, tiempo, 1]
                contador_frailty_especifico += 1
        if  contador_frailty_especifico == 0:
            curva_temporal.append(0)
        else:
            curva_temporal.append(valor_frailty_tiempo_especifico/contador_frailty_especifico)
    frailty_curves_general.append(curva_temporal)
frailty_curves_general = np.array(frailty_curves_general)
#%%
plt.scatter(np.arange(0,100),frailty_curves_general[5])
plt.scatter(np.arange(0,100),frailty_curves_general[-15])
#%%
#persona especifica, tiempo especifico, indicativo especifico
mortality_rate_modification[0][7, 90, 2]
# %%
frailty_index_curves = []
for celulas in mortality_rate_modification:

    suma = np.nansum(celulas[:, :, 1], axis=0)
    longitud_valida = np.sum(~np.isnan(celulas[:, :, 1]), axis=0)
    promedio_curva_frailty_index = np.divide(suma, longitud_valida, out=np.zeros_like(suma), where=longitud_valida != 0)
    frailty_index_curves.append(promedio_curva_frailty_index)
frailty_index_curves = np.array(frailty_index_curves)/100
# %%
datos = frailty_index_curves
valores_parametro = np.linspace(0, 1, len(frailty_index_curves))
cmap = cm.get_cmap('plasma')
fig, ax = plt.subplots(figsize = (8,5))

for i in range(len(frailty_index_curves)):
    color = cmap(valores_parametro[i])
    ax.plot(np.arange(100), datos[i], color=color, alpha=0.7, label=f'Curva {i + 1}')

norm = Normalize(vmin=0.01, vmax=0.9)  # Ajusta los límites del color según tus necesidades
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])  # Necesario para que funcione el colorbar
cbar = plt.colorbar(sm, ax=ax, label=r'Mortality Rate $\mu_0$')
# Etiquetas de los ejes
ax.set_title(r"Frailty Index vs time" + "\n" + r"with different mortality rate $\mu_0$", fontsize = 16)
ax.set_xlabel(r'Time', fontsize = 14)
ax.set_ylabel(r'Frailty Index', fontsize = 14)
plt.show()
#%%
mortality_data_mortality_rate = []
for celulas in mortality_rate_modification:
    mortality_data = []
    for i in np.arange(0,100):
        mortality_data.append(np.sum(celulas[:,i,2]))
    mortality_data = np.array(mortality_data)
    mortality_data_mortality_rate.append(mortality_data)
mortality_data_mortality_rate = np.array(mortality_data_mortality_rate)
#%%
datos = mortality_data_mortality_rate*3.3
valores_parametro = np.linspace(0, 1, len(mortality_data_mortality_rate))
cmap = cm.get_cmap('plasma')
fig, ax = plt.subplots(figsize = (8,5))

for i in range(len(mortality_data_mortality_rate)):
    color = cmap(valores_parametro[i])
    ax.plot(np.arange(100), datos[i], color=color, alpha=0.7, label=f'Curva {i + 1}')

norm = Normalize(vmin=0.01, vmax=0.9)  # Ajusta los límites del color según tus necesidades
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])  # Necesario para que funcione el colorbar
cbar = plt.colorbar(sm, ax=ax, label=r'Mortality Rate $\mu_0$')
# Etiquetas de los ejes
ax.set_title(r"Mortality $\mu$ vs time" + "\n" + r"with different mortality rate $\mu_0$", fontsize = 16)
ax.set_xlabel(r'Time', fontsize = 14)
ax.set_ylabel(r'Mortality $\mu$', fontsize = 14)
plt.show()
# %%
mortality_exponent_modification = np.load('/Users/rafther0112/Documents/GitHub/AGING_RESULTS_SIMULATIONS/Simulacion_modificacion_exponente_mortalidad.npy', allow_pickle=True)
# %%
frailty_index_curves_exponent = []
for celulas in (mortality_exponent_modification):

    suma = np.nansum(celulas[:, :, 1], axis=0)
    longitud_valida = np.sum(~np.isnan(celulas[:, :, 1]), axis=0)
    promedio_curva_frailty_index = np.divide(suma, longitud_valida, out=np.zeros_like(suma), where=longitud_valida != 0)
    frailty_index_curves_exponent.append(promedio_curva_frailty_index)
frailty_index_curves_exponent = np.array(frailty_index_curves_exponent)/100
#%%
datos = frailty_index_curves_exponent
valores_parametro = np.linspace(0, 1, len(frailty_index_curves_exponent))
cmap = cm.get_cmap('plasma')
fig, ax = plt.subplots(figsize=(8,5))
for i in range(len(frailty_index_curves_exponent)):
    color = cmap(valores_parametro[i])
    ax.plot(np.arange(100), datos[i], color=color, alpha=0.7, label=f'Curva {i + 1}')

norm = Normalize(vmin=1, vmax=10)  # Ajusta los límites del color según tus necesidades
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])  # Necesario para que funcione el colorbar
cbar = plt.colorbar(sm, ax=ax, label='Exponent of mortality C')

# Etiquetas de los ejes
ax.set_title(r"Frailty Index vs time" + "\n" + r"with different exponent of mortality  $C$", fontsize = 16)
ax.set_xlabel(r'Time', fontsize = 14)
ax.set_ylabel(r'Frailty Index', fontsize = 14)
plt.show()
# %%
mortality_data_exponent_mortality = []
for celulas in mortality_exponent_modification:
    mortality_data = []
    for i in np.arange(0,100):
        mortality_data.append(np.sum(celulas[:,i,2]))
    mortality_data = np.array(mortality_data)
    mortality_data_exponent_mortality.append(mortality_data)
mortality_data_exponent_mortality = np.array(mortality_data_exponent_mortality)
# %%
datos = mortality_data_exponent_mortality
valores_parametro = np.linspace(0, 1, len(frailty_index_curves_exponent))
cmap = cm.get_cmap('plasma')
fig, ax = plt.subplots(figsize=(8,5))
for i in range(len(frailty_index_curves_exponent)):
    color = cmap(valores_parametro[i])
    ax.plot(np.arange(100), datos[i], color=color, alpha=0.7, label=f'Curva {i + 1}')

norm = Normalize(vmin=1, vmax=10)  # Ajusta los límites del color según tus necesidades
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])  # Necesario para que funcione el colorbar
cbar = plt.colorbar(sm, ax=ax, label='Exponent of mortality C')

# Etiquetas de los ejes
ax.set_title(r"Mortality $\mu$ vs time" + "\n" + r"with different exponent of mortality  $C$", fontsize = 16)
ax.set_xlabel(r'Time', fontsize = 14)
ax.set_ylabel(r'Mortality $\mu$', fontsize = 14)
plt.show()
# %%
damage_rate_modification = np.load('/Users/rafther0112/Documents/GitHub/AGING_RESULTS_SIMULATIONS/Simulacion_modificacion_tasa_daño_exogeno_A.npy', allow_pickle=True)
#%%
frailty_index_curves_damage_rate = []
for celulas in damage_rate_modification:

    suma = np.nansum(celulas[:, :, 1], axis=0)
    longitud_valida = np.sum(~np.isnan(celulas[:, :, 1]), axis=0)
    promedio_curva_frailty_index = np.divide(suma, longitud_valida, out=np.zeros_like(suma), where=longitud_valida != 0)
    frailty_index_curves_damage_rate.append(promedio_curva_frailty_index)
frailty_index_curves_damage_rate = np.array(frailty_index_curves_damage_rate)/100
#%%
datos = frailty_index_curves_damage_rate
valores_parametro = np.linspace(0, 1, len(frailty_index_curves_damage_rate))
cmap = cm.get_cmap('plasma')
fig, ax = plt.subplots(figsize = (8,5))

for i in range(len(frailty_index_curves_damage_rate)):
    color = cmap(valores_parametro[i])
    ax.plot(np.arange(100), datos[i], color=color, alpha=0.7, label=f'Curva {i + 1}')

norm = Normalize(vmin=0.001, vmax=0.1)  # Ajusta los límites del color según tus necesidades
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])  # Necesario para que funcione el colorbar
cbar = plt.colorbar(sm, ax=ax, label=r'Damage rate $a$')
# Etiquetas de los ejes
ax.set_title(r"Frailty Index vs time" + "\n" + r"with different damage rate $a$", fontsize = 16)
ax.set_xlabel(r'Time', fontsize = 14)
ax.set_ylabel(r'Frailty Index', fontsize = 14)
plt.show()
#%%
mortality_data_damage_rate = []
for celulas in damage_rate_modification:
    mortality_data = []
    for i in np.arange(0,100):
        mortality_data.append(np.sum(celulas[:,i,2]))
    mortality_data = np.array(mortality_data)
    mortality_data_damage_rate.append(mortality_data)
mortality_data_damage_rate = np.array(mortality_data_damage_rate)
#%%
datos = mortality_data_damage_rate
valores_parametro = np.linspace(0, 1, len(mortality_data_damage_rate))
cmap = cm.get_cmap('plasma')
fig, ax = plt.subplots(figsize = (8,5))

for i in range(len(mortality_data_damage_rate)):
    color = cmap(valores_parametro[i])
    ax.plot(np.arange(100), datos[i], color=color, alpha=0.7, label=f'Curva {i + 1}')

norm = Normalize(vmin=0.001, vmax=0.1)  # Ajusta los límites del color según tus necesidades
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])  # Necesario para que funcione el colorbar
cbar = plt.colorbar(sm, ax=ax, label=r'Damage rate $a$')
# Etiquetas de los ejes
ax.set_title(r"Mortality vs time" + "\n" + r"with different damage rate $a$", fontsize = 16)
ax.set_xlabel(r'Time', fontsize = 14)
ax.set_ylabel(r'Frailty Index', fontsize = 14)
plt.show()
# %%
recovery_rate_modification = np.load('/Users/rafther0112/Documents/GitHub/AGING_RESULTS_SIMULATIONS/Simulacion_modificacion_tasa_reparación_R_nuevo.npy', allow_pickle=True)
#%%
frailty_index_curves_recovery_rate = []
for celulas in recovery_rate_modification:

    suma = np.nansum(celulas[:, :, 1], axis=0)
    longitud_valida = np.sum(~np.isnan(celulas[:, :, 1]), axis=0)
    promedio_curva_frailty_index = np.divide(suma, longitud_valida, out=np.zeros_like(suma), where=longitud_valida != 0)
    frailty_index_curves_recovery_rate.append(promedio_curva_frailty_index)
frailty_index_curves_recovery_rate = np.array(frailty_index_curves_recovery_rate)/100
#%%
datos = frailty_index_curves_recovery_rate
valores_parametro = np.linspace(0, 1, len(frailty_index_curves_recovery_rate))
cmap = cm.get_cmap('plasma')
fig, ax = plt.subplots(figsize = (8,5))

for i in range(len(frailty_index_curves_recovery_rate)):
    color = cmap(valores_parametro[i])
    ax.plot(np.arange(100), datos[i], color=color, alpha=0.7, label=f'Curva {i + 1}')

norm = Normalize(vmin=0.1, vmax=1)  # Ajusta los límites del color según tus necesidades
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])  # Necesario para que funcione el colorbar
cbar = plt.colorbar(sm, ax=ax, label=r'Recovery rate $r$')
# Etiquetas de los ejes
ax.set_title(r"Frailty Index vs time" + "\n" + r"with different recovery rate $r$", fontsize = 16)
ax.set_xlabel(r'Time', fontsize = 14)
ax.set_ylabel(r'Frailty Index', fontsize = 14)
plt.show()
#%%
mortality_data_recovery_rate = []
for celulas in recovery_rate_modification:
    mortality_data = []
    for i in np.arange(0,100):
        mortality_data.append(np.sum(celulas[:,i,2]))
    mortality_data = np.array(mortality_data)
    mortality_data_recovery_rate.append(mortality_data)
mortality_data_recovery_rate = np.array(mortality_data_recovery_rate)
# %%
datos = mortality_data_recovery_rate
valores_parametro = np.linspace(0, 1, len(mortality_data_recovery_rate))
cmap = cm.get_cmap('plasma')
fig, ax = plt.subplots(figsize = (8,5))

for i in range(len(mortality_data_recovery_rate)):
    color = cmap(valores_parametro[i])
    ax.plot(np.arange(100), datos[i], color=color, alpha=0.7, label=f'Curva {i + 1}')

norm = Normalize(vmin=0.1, vmax=1)  # Ajusta los límites del color según tus necesidades
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])  # Necesario para que funcione el colorbar
cbar = plt.colorbar(sm, ax=ax, label=r'Recovery rate $r$')
# Etiquetas de los ejes
ax.set_title(r"Mortality vs time" + "\n" + r"with different recovery rate $r$", fontsize = 16)
ax.set_xlabel(r'Time', fontsize = 14)
ax.set_ylabel(r'Frailty Index', fontsize = 14)
plt.show()
# %%
number_nodes_modification = np.load('/Users/rafther0112/Documents/GitHub/AGING_RESULTS_SIMULATIONS/Simulacion_modificacion_numero_nodos.npy', allow_pickle=True)
#%%
frailty_index_number_nodes_modification= []
for celulas in number_nodes_modification:

    suma = np.nansum(celulas[:, :, 1], axis=0)
    longitud_valida = np.sum(~np.isnan(celulas[:, :, 1]), axis=0)
    promedio_curva_frailty_index = np.divide(suma, longitud_valida, out=np.zeros_like(suma), where=longitud_valida != 0)
    frailty_index_number_nodes_modification.append(promedio_curva_frailty_index)
frailty_index_number_nodes_modification = np.array(frailty_index_number_nodes_modification)
# %%
datos = frailty_index_number_nodes_modification
valores_parametro = np.linspace(0, 1, len(frailty_index_number_nodes_modification))
cmap = cm.get_cmap('plasma')
fig, ax = plt.subplots(figsize = (8,5))
valores_de_nodos = [1,5,10,15,20,25,50,100,1000]
for i in range(len(frailty_index_number_nodes_modification)):
    color = cmap(valores_parametro[i])
    ax.plot(np.arange(100), datos[i]/valores_de_nodos[i], color=color, alpha=0.7, label=f'Curva {i + 1}')

norm = Normalize(vmin=1, vmax=1000)  # Ajusta los límites del color según tus necesidades
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])  # Necesario para que funcione el colorbar
cbar = plt.colorbar(sm, ax=ax, label=r'Damage rate $a$')
# Etiquetas de los ejes
ax.set_title(r"Frailty Index vs time" + "\n" + r"with different number of nodes $N_{t}$", fontsize = 16)
ax.set_xlabel(r'Time', fontsize = 14)
ax.set_ylabel(r'Frailty Index', fontsize = 14)
plt.show()
# %%
mortality_data_number_nodes= []
for celulas in number_nodes_modification:
    mortality_data = []
    for i in np.arange(0,100):
        mortality_data.append(np.sum(celulas[:,i,2]))
    mortality_data = np.array(mortality_data)
    mortality_data_number_nodes.append(mortality_data)
mortality_data_number_nodes = np.array(mortality_data_recovery_rate)
# %%
datos = mortality_data_number_nodes
valores_parametro = np.linspace(0, 1, len(mortality_data_number_nodes))
cmap = cm.get_cmap('plasma')
fig, ax = plt.subplots(figsize = (8,5))

for i in range(len(mortality_data_number_nodes)):
    color = cmap(valores_parametro[i])
    ax.plot(np.arange(100), datos[i], color=color, alpha=0.7, label=f'Curva {i + 1}')

norm = Normalize(vmin=1, vmax=1000)  # Ajusta los límites del color según tus necesidades
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])  # Necesario para que funcione el colorbar
cbar = plt.colorbar(sm, ax=ax, label=r'Damage rate $a$')
# Etiquetas de los ejes
ax.set_title(r"Mortality vs time" + "\n" + r"with different number of nodes rate $N$", fontsize = 16)
ax.set_xlabel(r'Time', fontsize = 14)
ax.set_ylabel(r'Frailty Index', fontsize = 14)
plt.show()
# %%
