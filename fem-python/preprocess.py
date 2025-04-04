import numpy as np

class Mesh:
    def __init__(self, x, y, z, IEN, ID, LM, ned, neq, nnp, ng, BC_g_list, BC_fix_list, free_range, freefix_range):
        self.x = x
        self.y = y
        self.z = z
        self.IEN = IEN
        self.ID = ID
        self.LM = LM
        self.ned = ned
        self.neq = neq
        self.nnp = nnp
        self.ng = ng
        self.BC_g_list = BC_g_list
        self.BC_fix_list = BC_fix_list
        self.free_range = free_range
        self.freefix_range = freefix_range

class ElementDef:
    def __init__(self, name, nen, internal_dim, gmsh_number):
        self.name = name
        self.nen = nen
        self.internal_dim = internal_dim
        self.gmsh_number = gmsh_number

# Standard GMSH element definitions
std_element_defs = {
    "line": ElementDef("Linear Line", 2, 1, 1),
    "triangle": ElementDef("Linear Triangle", 3, 2, 2),
    "quad": ElementDef("4 node quadrilateral", 4, 2, 3)
}

def build_ID(nnp, g_list, ned, fix_list):
    neq = 0
    count = 0
    
    # Compute ng
    ng = np.sum(fix_list > 0)
    
    # Construct ID
    ID = np.zeros((ned, nnp), dtype=int)
    totaldof = ned * nnp
    gg = np.zeros(ng) if ng > 0 else []
    
    if ng > 0:
        count = 0
        for A in range(nnp):
            for i in range(ned):
                if fix_list[i, A] == 1:
                    ID[i, A] = totaldof - count
                    gg[count] = g_list[i, A]
                    count += 1
                else:
                    neq += 1
                    ID[i, A] = neq
    else:
        for A in range(nnp):
            for i in range(ned):
                neq += 1
                ID[i, A] = neq
        gg = []

    # Fix numbering to be zero-based
    ID = ID - 1
    
    free_range = range(neq)
    freefix_range = range(neq, neq + ng)
    
    return ID, neq, ng, free_range, freefix_range

def build_LM(ID, IEN):
    ned = ID.shape[0]
    LM = {}
    
    for element_type, ien in IEN.items():
        nel = ien.shape[0]
        nen = std_element_defs[element_type].nen
        lm = np.zeros((nen * ned, nel), dtype=int)
        
        for e in range(nel):
            for a in range(nen):
                for i in range(ned):
                    p = a + nen * (i)
                    lm[p, e] = ID[i, ien[e, a]]
        
        LM[element_type] = lm
    
    return LM

def build_mesh(x, y, z, IEN, DofPerNode, BC_fix_list, BC_g_list):
    nnp = len(x)
    ID, neq, ng, free_range, freefix_range = build_ID(nnp, BC_g_list, DofPerNode, BC_fix_list)
    LM = build_LM(ID, IEN)
    
    return Mesh(x, y, z, IEN, ID, LM, DofPerNode, neq, nnp, ng, BC_g_list, BC_fix_list, free_range, freefix_range)
