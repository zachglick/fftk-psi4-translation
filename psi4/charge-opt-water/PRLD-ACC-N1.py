import psi4
import optking
import qcelemental as qcel
import numpy as np

# Summary: optimizes intermolecular coordinates between two frozen monomers.

# Split dimers in input with "--"
mol = psi4.geometry("""
    N1 0.843999981880188 0.4560000002384186 0.9229999780654907
    H2 1.371999979019165 0.4950000047683716 1.7979999780654907
    C3 1.6430000066757202 1.125 -0.11299999803304672
    H4 2.0 2.0969998836517334 0.24300000071525574
    H5 0.9869999885559082 1.3070000410079956 -0.9710000157356262
    C6 2.7790000438690186 0.15299999713897705 -0.5040000081062317
    H7 3.688999891281128 0.3840000033378601 0.05900000035762787
    H8 3.0309998989105225 0.22200000286102295 -1.565999984741211
    C9 2.2309999465942383 -1.2400000095367432 -0.1120000034570694
    H10 2.865999937057495 -1.7070000171661377 0.6470000147819519
    H11 2.184999942779541 -1.9279999732971191 -0.9620000123977661
    C12 0.8299999833106995 -0.9409999847412109 0.4690000116825104
    H13 0.5419999957084656 -1.6089999675750732 1.2860000133514404
    H14 0.06800000369548798 -1.0269999504089355 -0.3140000104904175
    --
    H       1     rAH       2  112.65       3  -124.88
    x      15     1.0       1   90.00       2    0.00
    O      15  0.9572       16   90.00       1  180.00
    H      17  0.9572      15  104.52        16     dih

    rAH  = 2.0
    dih = 0.0 
    nocom
    unit angstrom
    """)

# Ordering of reference atoms is essential for properly defining contraints.
# Note that if either of the theta angles are linear. optking will crash.
# Keep to the ordering below and everything should be fine :)
#  A1: Ligand donor/acceptor atom
#  A2: Ligand atom in-plane with dih (NOT in-line with H-bond)
#  A3: Ligand atom
#  B1: Water donor/acceptor atom
#  B2: Water atom in-plane with dih (NOT in-line with H-bond)
#  B3: Final water atom

# Here are the dimer coordinates that are used with their definitions.
#  R         : Distance A1 to B1 (rAH)
#  theta_A   : Angle,          A2-A1-B1
#  theta_B   : Angle,          A1-B1-B2
#  tau       : Dihedral angle, A2-A1-B1-B2 (dih)
#  phi_A     : Dihedral angle, A3-A2-A1-B1
#  phi_B     : Dihedral angle, A1-B1-B2-B3

# Coords are 1-indexed without dummies here
dimer = {
    "Natoms per frag": [14, 3],
    "A Frag": 1,
    "A Ref Atoms": [[1],[2],[3]],
    "A Label": "PRLD",
    "B Frag": 2,
    "B Ref Atoms": [[15],[17],[16]],
    "B Label": "Water",
    "Frozen": ["theta_A", "theta_B", "phi_A", "phi_B"],
}

psi4.set_memory("1.0 GB")
psi4.set_output_file('PRLD-ACC-H1.out', False)
psi4.core.clean_options()
psi4_options = {
    "basis": "6-31G*",
    "d_convergence": 9,
    "frag_mode": "multi",
    "freeze_intrafrag": 'true',
}
psi4.set_options(psi4_options)
result = optking.optimize_psi4("hf", **{"interfrag_coords": str(dimer)})

# Can use qcel to output final rAH and dih values
xyzs = result["final_molecule"]["geometry"]
xyzs = np.array(xyzs)
xyzs = np.reshape(xyzs, (-1,3))
# coords are zero-indexed here:
rAH = qcel.util.measure_coordinates(xyzs,[0,14],True) # in bohr
dih = qcel.util.measure_coordinates(xyzs,[1,0,14,16],True) # in degrees

bohr_to_ang = qcel.constants.conversion_factor("bohr", "angstrom")
print(rAH * bohr_to_ang)
print(dih)
