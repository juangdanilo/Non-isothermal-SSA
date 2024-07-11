from ssa_simulation import SSA

params = {
    'k1': 1e-4, 'k2': 8e-6, 'k3': 1e-11, 'k4': 6e-5, 'k5': 3e-4, 'k6': 4e-4,
    'A0': 500, 'B0': 500, 'C0': 99000, 'Cj0': 0, 'P0': 0, 'T0': 50,
    'E1': 353, 'E2': 353, 'E3': 353, 'E4': 353, 'E5': 353, 'E6': 353,
    'beta': 10, 'q': 1.5, 'Nc': 100, 'shape': 1000
}

ssa = SSA(params)

result = ssa.organize()

print(result)