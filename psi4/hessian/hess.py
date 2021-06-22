import numpy as np
import psi4
import optking

# Edit these
psi4.core.set_output_file('hess.out')
psi4.set_memory("1GB")
psi4.set_num_threads(6)

# Specify geometry here
mol = psi4.geometry("""
0 1
N           -0.851383575847     0.545401576222     0.747813415726
H           -0.352345727928     0.599592947948     1.638193413524
C           -0.026593266302     1.207924140237    -0.269928358509
H            0.308577012235     2.187043048549     0.082675199158
H           -0.656712212286     1.370738247192    -1.151465407112
C            1.128794127761     0.242758496768    -0.610932136414
H            2.016952054235     0.488920787753    -0.021166045504
H            1.415796771435     0.300401418793    -1.664791753657
C            0.580896761378    -1.149189703140    -0.219935511754
H            1.195126521620    -1.600607422633     0.565033996747
H            0.569407550630    -1.847642807610    -1.061566673206
C           -0.838702968673    -0.855996638817     0.310812574256
H           -1.148362789281    -1.516734500287     1.124984603326
H           -1.573080054265    -0.957331636962    -0.496081473508
""")

# Specify internal coordinates  here
# Note: atoms indices start at zero
bonds = [
(1, 0),
(2, 0),
(11, 0),
(3, 2),
(4, 2),
(5, 2),
(6, 5),
(7, 5),
(8, 5),
(9, 8),
(10, 8),
(11, 8),
(12, 11),
(13, 11),
]

angles = [
(3, 2, 0),
(4, 2, 0),
(5, 2, 0),
(8, 11, 0),
(12, 11, 0),
(13, 11, 0),
(2, 0, 1),
(11, 0, 1),
(11, 0, 2),
(6, 5, 2),
(7, 5, 2),
(8, 5, 2),
(4, 2, 3),
(5, 2, 3),
(5, 2, 4),
(9, 8, 5),
(10, 8, 5),
(11, 8, 5),
(7, 5, 6),
(8, 5, 6),
(8, 5, 7),
(12, 11, 8),
(13, 11, 8),
(10, 8, 9),
(11, 8, 9),
(11, 8, 10),
(13, 11, 12),
]

dihedrals = [
(3, 2, 0, 1),
(4, 2, 0, 1),
(5, 2, 0, 1),
(3, 2, 0, 11),
(4, 2, 0, 11),
(5, 2, 0, 11),
(8, 11, 0, 1),
(12, 11, 0, 1),
(13, 11, 0, 1),
(8, 11, 0, 2),
(12, 11, 0, 2),
(13, 11, 0, 2),
(6, 5, 2, 0),
(7, 5, 2, 0),
(8, 5, 2, 0),
(6, 5, 2, 3),
(7, 5, 2, 3),
(8, 5, 2, 3),
(6, 5, 2, 4),
(7, 5, 2, 4),
(8, 5, 2, 4),
(9, 8, 5, 2),
(10, 8, 5, 2),
(11, 8, 5, 2),
(9, 8, 5, 6),
(10, 8, 5, 6),
(11, 8, 5, 6),
(9, 8, 5, 7),
(10, 8, 5, 7),
(11, 8, 5, 7),
(0, 11, 8, 5),
(12, 11, 8, 5),
(13, 11, 8, 5),
(0, 11, 8, 9),
(12, 11, 8, 9),
(13, 11, 8, 9),
(0, 11, 8, 10),
(12, 11, 8, 10),
(13, 11, 8, 10),
]

psi4_options  = {
"mp2_type" : "df",
"freeze_core" : True,
"basis" : "6-31G*",
"g_convergence" : "gau_tight",
}

psi4.set_options(psi4_options)

# Optimize and extract information, inc. final geometry
json_output = optking.optimize_psi4('MP2')
E   = json_output['energies'][-1]
print(f"Optimized Energy: {E}")
xyz_array = np.array(json_output['final_molecule']['geometry'])
xyz = xyz_array.reshape(mol.natom(),3)
mol.set_geometry(psi4.core.Matrix.from_array(xyz))

optking.optparams.Params = optking.optparams.OptParams({})

from optking import stre, bend, tors
bonds = [stre.Stre(*bond) for bond in bonds]
angles = [bend.Bend(*angle) for angle in angles]
dihedrals = [tors.Tors(*dihedral) for dihedral in dihedrals]
coords = bonds + angles + dihedrals

psi4.set_options(psi4_options)
Z = [mol.Z(i) for i in range(0,mol.natom())]
masses = [mol.mass(i) for i in range(0,mol.natom())]
f1 = optking.frag.Frag(Z, xyz, masses, intcos=coords, frozen=False)
OptMol = optking.molsys.Molsys([f1])
print(OptMol)

# Compute the Cartesian Hessian with psi4
Hxyz = psi4.hessian('MP2') # returns a psi4 matrix
print("Cartesian hessian")
print(Hxyz.to_array())

# Transform hessian to internals with optking, returns a numpy array
Hq = optking.intcosMisc.hessian_to_internals(Hxyz.to_array(), OptMol)
print("Internal Coordinate Hessian")
print(Hq)

# Also print transformed Hessian to output File
psi4.core.Matrix.from_array(Hq).print_out()

# Note: there are probably better ways to print out the internal coordinate hessian that would minimize the amount of parsing needed in FFTK
