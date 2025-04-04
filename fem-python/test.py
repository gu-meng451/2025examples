#%%
import preprocess as pp
import numpy as np

# %%
# make an element
nnp = 5
nel = 4
nee = 2
# make dictionary of elements
IEN = {"line": np.zeros((nel, nee), dtype=int)}
IEN["line"][0, :] = np.array([1, 3]) - 1
IEN["line"][1, :] = np.array([3, 4]) - 1
IEN["line"][2, :] = np.array([4, 5]) - 1
IEN["line"][3, :] = np.array([5, 2]) - 1

# %%
x = np.linspace(0, 1, nnp)
x = x[[0,4,1,2,3]]

# %%
BC_fix_list = np.zeros((1,nnp), dtype=int)
BC_fix_list[0, 0] = 1

BC_g_list = np.zeros((1,nnp), dtype=int)
BC_g_list[0, 0] = 1

# %%
mesh = pp.build_mesh(x, [], [], IEN, 1, BC_fix_list, BC_g_list)
# %%
