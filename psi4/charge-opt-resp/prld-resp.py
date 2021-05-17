import psi4
import resp

mol = psi4.geometry("""
    0 1
    N       0.84399998       0.45600000       0.92299998
    H       1.37199998       0.49500000       1.79799998
    C       1.64300001       1.12500000      -0.11300000
    H       2.00000000       2.09699988       0.24300000
    H       0.98699999       1.30700004      -0.97100002
    C       2.77900004       0.15300000      -0.50400001
    H       3.68899989       0.38400000       0.05900000
    H       3.03099990       0.22200000      -1.56599998
    C       2.23099995      -1.24000001      -0.11200000
    H       2.86599994      -1.70700002       0.64700001
    H       2.18499994      -1.92799997      -0.96200001
    C       0.82999998      -0.94099998       0.46900001
    H       0.54200000      -1.60899997       1.28600001
    H       0.06800000      -1.02699995      -0.31400001
    """)
mol.update_geometry()

options = { 'RESP_A' : 0.0005 }
charges_1 = resp.resp([mol], [options])
print(charges_1[0][1])
resp.set_stage2_constraint(mol, charges_1[0][1], options)
charges_2 = resp.resp([mol], [options])

# Get RESP charges
print("\nStage Two:\n")
print('RESP Charges')
print(charges_2[0][1])
