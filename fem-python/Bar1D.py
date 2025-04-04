import numpy as np
from shapefunctions import shapefunc
from preprocess import std_element_defs  # Assuming this exists in Preprocess

def assemble_stiffness(mesh, properties):
    """Assemble the global stiffness matrix."""
    ned = mesh.ned
    totaldofs = ned * mesh.nnp
    K = np.zeros((totaldofs, totaldofs))

    # Loop over each element in the mesh
    for element_type, ien in mesh.IEN.items():
        nel = ien.shape[0]
        nen = std_element_defs[element_type].nen
        nee = ned * nen
        LM = mesh.LM[element_type]

        for e in range(nel):
            A = ien[e, :nen]
            xe = mesh.x[A]

            ke = element_stiffness(xe, properties)

            # Assemble element stiffness into the global stiffness matrix
            for loop1 in range(nee):
                i = LM[loop1, e]
                for loop2 in range(nee):
                    j = LM[loop2, e]
                    K[i, j] += ke[loop1, loop2]

    return K

def element_stiffness(xe, properties):
    """Compute the element stiffness matrix for a 1D bar."""
    Le = abs(xe[1] - xe[0])
    E = properties.E
    A = properties.A

    ke = (E * A / Le) * np.array([[1.0, -1.0],
                                   [-1.0, 1.0]])

    return ke

def assemble_rhs(mesh, external_forcing, quad_rules):
    """Assemble the global right-hand-side force vector."""
    ned = mesh.ned
    totaldofs = ned * mesh.nnp
    F = np.zeros(totaldofs)

    # Loop over each element in the mesh
    for element_type, ien in mesh.IEN.items():
        nel = ien.shape[0]
        N = shapefunc(element_type)
        element_quad_rule = quad_rules[element_type]
        nen = std_element_defs[element_type].nen
        nee = ned * nen
        LM = mesh.LM[element_type]

        for e in range(nel):
            A = ien[e, :nen]
            xe = mesh.x[A]

            fe = element_forcing(xe, N, element_quad_rule, external_forcing)

            # Assemble element force into the global force vector
            for loop1 in range(nee):
                i = LM[loop1, e]
                F[i] += fe[loop1]

    return F

def element_forcing(xe, N, element_quad_rule, external_forcing):
    """Compute the element force vector."""
    ned = 1
    nen = len(xe)
    nee = ned * nen
    fe = np.zeros(nee)

    # Integration loop
    for ξ, w in element_quad_rule.iterator:
        # Evaluate the shape function
        Ne, Nξ = N(ξ)

        # Evaluate the external loading at x(ξ)
        x = np.dot(Ne, xe)
        fext = external_forcing(x)

        # Build Jacobian
        detJ = np.dot(Nξ, xe)

        # Compute dV0
        dV0 = detJ * w

        # Integrate
        fe += Ne * fext * dV0

    return fe
