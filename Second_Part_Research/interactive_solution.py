#%%
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import ipywidgets as widgets
from IPython.display import display

def frailty_index_differential_equation(f, t, a, b, r, s):
    dfdt = a * (1 - f) * (1 + b * t) - f * r * (1 - s * t)
    return dfdt

# Initial conditions
f0 = 0.5  # Initial frailty index
a = 1.2
b = 0.5
r = 0.8
s = 0.2

# Time points
t = np.linspace(0, 10, 100)  # 100 time steps from 0 to 10

# Create figure and axis for the plot
fig, ax = plt.subplots(figsize=(8, 6))

# Function to update the plot
def update_plot(a, b, r, s):
    f_solution = odeint(frailty_index_differential_equation, f0, t, args=(a, b, r, s))
    
    ax.clear()  # Clear the previous plot
    ax.plot(t, f_solution, label='Frailty Index')
    ax.set_xlabel('Time')
    ax.set_ylabel('Frailty Index')
    ax.set_title('Solution of Frailty Index Differential Equation')
    ax.legend()
    ax.grid(True)
    plt.show()

# Create sliders for interactive parameter adjustment
slider_a = widgets.FloatSlider(value=a, min=0.1, max=2.0, step=0.1, description='a:')
slider_b = widgets.FloatSlider(value=b, min=0, max=2, step=0.1, description='b:')
slider_r = widgets.FloatSlider(value=r, min=0, max=1, step=0.1, description='r:')
slider_s = widgets.FloatSlider(value=s, min=0, max=1, step=0.1, description='s:')

# Define the function to be called when sliders are changed
def on_value_change(change):
    a = slider_a.value
    b = slider_b.value
    r = slider_r.value
    s = slider_s.value
    update_plot(a, b, r, s)

# Attach the function to the slider value change event
slider_a.observe(on_value_change, names='value')
slider_b.observe(on_value_change, names='value')
slider_r.observe(on_value_change, names='value')
slider_s.observe(on_value_change, names='value')

# Display the interactive plot and sliders
display(widgets.VBox([slider_a, slider_b, slider_r, slider_s]))
update_plot(a, b, r, s)  # Display the initial plot

# %%
