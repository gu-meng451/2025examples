import numpy as np
import matplotlib.pyplot as plt
from shapefunctions import shapefunc

def init_plot(xlabel="Position x", ylabel="Displacement u"):
    """Initialize a plot with given labels."""
    fig, ax = plt.subplots()
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    return fig, ax

def draw_element(ax, mesh, q, linestyle='-', linewidth=2, color='C0'):
    """Draw elements with their displacements."""
    for element_type, ien in mesh.IEN.items():
        if element_type == "line":
            nel = ien.shape[0]
            for e in range(nel):
                A = ien[e, :]
                xpts = mesh.x[A]
                idx = mesh.ID[0, A]  # Adjusting for Python's 0-based indexing
                qe = q[idx]
                ax.plot(xpts, qe, linestyle=linestyle, linewidth=linewidth, color=color)

def draw_element_stress(ax, mesh, q, properties, linestyle='-', linewidth=2, color='C0'):
    """Draw elements with stress distribution."""
    for element_type, ien in mesh.IEN.items():
        if element_type == "line":
            nel = ien.shape[0]
            N = shapefunc(element_type)
            ξ_vals = np.linspace(-1, 1, 3)

            for e in range(nel):
                A = ien[e, :]
                xe = mesh.x[A]
                idx = mesh.ID[0, A]
                qe = q[idx]

                X = np.zeros(len(ξ_vals))
                σ = np.zeros(len(ξ_vals))

                for i, ξ in enumerate(ξ_vals):
                    NN, Nξ = N(ξ)
                    X[i] = np.dot(NN, xe)
                    σ[i] = properties.E * np.dot(Nξ, qe) / np.dot(Nξ, xe)

                ax.plot(X, σ, linestyle=linestyle, linewidth=linewidth, color=color)
