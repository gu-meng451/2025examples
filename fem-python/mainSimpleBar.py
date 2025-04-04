# %%
import numpy as np
from collections import namedtuple
from preprocess import build_mesh
import quadrature as quad
from Bar1D import assemble_stiffness, assemble_rhs
import simpleVisualization as vis
import matplotlib.pyplot as plt

# Define properties
Properties = namedtuple("Properties", ["E", "A", "L", "f0"])
prop = Properties(E=1.0, A=1.0, L=1.0, f0=1.0)

# Define quad rule
quad_rules = {"line": quad.gauss_legendre_1d(2)}

# Define the mesh
nnp = 5
nel = 4
nee = 2
# make dictionary of elements
IEN = {"line": np.zeros((nel, nee), dtype=int)}
IEN["line"][0, :] = np.array([1, 3]) - 1
IEN["line"][1, :] = np.array([3, 4]) - 1
IEN["line"][2, :] = np.array([4, 5]) - 1
IEN["line"][3, :] = np.array([5, 2]) - 1

# rearrange to match the class-ordered problem
x = np.linspace(0, 1, nnp)
x = x[[0, 4, 1, 2, 3]]

# flag the essential boundary conditions
BC_fix_list = np.zeros((1, nnp), dtype=int)
BC_fix_list[0, 0] = 1

# values for the essential boundary conditions
BC_g_list = np.zeros((1, nnp), dtype=int)
BC_g_list[0, 0] = 0.0

mesh = build_mesh(x, [], [], IEN, 1, BC_fix_list, BC_g_list)


# Loading force function
def f(x):
    return prop.f0


# %% Assemble the global stiffness matrix
K = assemble_stiffness(mesh, prop)

# %% Assemble the global right-hand-side force vector
F = assemble_rhs(mesh, f, quad_rules)

# %% Solve
q = np.zeros(mesh.nnp * mesh.ned)
r1 = mesh.free_range
r2 = mesh.freefix_range
q[r1] = np.linalg.solve(K[np.ix_(r1, r1)], F[r1] - K[np.ix_(r1, r2)] @ q[r2])

# %%
fix, ax = vis.init_plot(xlabel="Position x", ylabel="Displacement u")
vis.draw_element(ax, mesh, q, linestyle="--", linewidth=2, color="r")


def u_true(x, prop):
    return prop.f0 / (2 * prop.E * prop.A) * (2 * prop.L - x) * x


x_vals = np.linspace(0, prop.L, 20)
u_vals = u_true(x_vals, prop)
ax.plot(
    x_vals, u_vals, label="Analytical Solution", color="C0", linewidth=2, linestyle="-"
)
# Show the figure
plt.legend()
plt.show()


# %%
def σ_true(x, prop):
    return prop.f0 / prop.A * (prop.L - x)


fig, ax = vis.init_plot(ylabel="Stress")
vis.draw_element_stress(ax, mesh, q, prop, linestyle="--", linewidth=2, color="r")
x_vals = np.linspace(0, prop.L, 20)
σ_vals = σ_true(x_vals, prop)
ax.plot(
    x_vals, σ_vals, label="Analytical Solution", color="C0", linewidth=2, linestyle="-"
)
# Show the figure
plt.legend()
plt.show()
# %%
