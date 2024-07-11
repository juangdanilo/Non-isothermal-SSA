import numpy as np
import pandas as pd
from math import exp

class SSA:
    def __init__(self, params):
        self.params = params
        self.initialize_arrays()

    def initialize_arrays(self):
        self.shape = self.params['shape']
        self.A = np.empty(self.shape + 1, dtype=float)
        self.B = np.empty(self.shape + 1, dtype=float)
        self.C = np.empty(self.shape + 1, dtype=float)
        self.Aj = np.empty(self.shape + 1, dtype=float)
        self.Bj = np.empty(self.shape + 1, dtype=float)
        self.Cj = np.empty(self.shape + 1, dtype=float)
        self.P = np.empty(self.shape + 1, dtype=float)
        self.D = np.empty(self.shape + 1, dtype=float)
        self.time = np.empty(self.shape + 1, dtype=float)
        self.T = np.empty(self.shape + 1, dtype=float)
        self.dif = np.empty(self.shape + 1, dtype=float)
        self.sum_propen = np.empty(self.shape + 1, dtype=float)
        self.propen_1 = np.empty(self.shape + 1, dtype=float)
        self.propen_2 = np.empty(self.shape + 1, dtype=float)
        self.propen_3 = np.empty(self.shape + 1, dtype=float)
        self.propen_4 = np.empty(self.shape + 1, dtype=float)
        self.propen_5 = np.empty(self.shape + 1, dtype=float)
        self.propen_6 = np.empty(self.shape + 1, dtype=float)

        self.A[0] = self.params['A0']
        self.B[0] = self.params['B0']
        self.C[0] = self.params['C0']
        self.Cj[0] = self.params['Cj0']
        self.P[0] = self.params['P0']
        self.T[0] = self.params['T0']
        self.Aj[0] = 0
        self.Bj[0] = 0
        self.D[0] = 0
        self.time[0] = 0
        self.dif[0] = 0
        self.sum_propen[0] = 0

    def run_simulation(self):
        for i in range(self.shape):
            self.step(i)

    def step(self, i):
        F = self.arrhenius(i)
        r1, r2 = np.random.random(2)
        propensities, sum_propen = self.calculate_propensities(F, i)
        self.update_state(r1, r2, propensities, sum_propen, i)

    def arrhenius(self, i):
        T_i = self.T[i]
        return [
            self.params['k1'] * exp(-self.params['E1'] / T_i),
            self.params['k2'] * exp(-self.params['E2'] / T_i),
            self.params['k3'] * exp(-self.params['E3'] / T_i),
            self.params['k4'] * exp(-self.params['E4'] / T_i),
            self.params['k5'] * exp(-self.params['E5'] / T_i),
            self.params['k6'] * exp(-self.params['E6'] / T_i)
        ]

    def calculate_propensities(self, F, i):
        q = self.params['q']
        A_i, B_i, C_i, Aj_i, Bj_i, Cj_i, P_i = self.A[i], self.B[i], self.C[i], self.Aj[i], self.Bj[i], self.Cj[i], self.P[i]
        propensities = np.array([
            F[0] * (A_i ** q) * (B_i ** q) * (self.params['A0'] ** (1 - q)) * (self.params['B0'] ** (1 - q)),
            F[1] * (Aj_i ** q) * (C_i ** q) * (self.params['C0'] ** (1 - q)),
            F[2] * (Bj_i ** q) * (C_i ** (2 * q)) * (self.params['C0'] ** (1 - q)),
            F[3] * (P_i ** q) * (B_i ** q) * (self.params['B0'] ** (1 - q)),
            F[4] * (Cj_i ** q) * (A_i ** q) * (self.params['A0'] ** (1 - q)),
            F[5] * (P_i ** q) * (Cj_i ** q)
        ])
        return propensities, propensities.sum()

    def update_state(self, r1, r2, propensities, sum_propen, i):
        if sum_propen > 0:
            r_new = (1 - (2 - self.params['q']) * r1) ** ((1 - self.params['q']) / (2 - self.params['q']))
            tau = (1 - r_new) / ((1 - self.params['q']) * sum_propen)

            cumulative_propensities = np.cumsum(propensities) / sum_propen
            reaction = np.searchsorted(cumulative_propensities, r2)

            self.update_particles(reaction, i)
            self.time[i + 1] = self.time[i] + tau
            self.T[i + 1] = self.params['T0'] + self.params['beta'] * self.time[i + 1]
            self.dif[i + 1] = 1 / tau
            self.sum_propen[i + 1] = sum_propen
            self.propen_1[i + 1] = propensities[0]
            self.propen_2[i + 1] = propensities[1]
            self.propen_3[i + 1] = propensities[2]
            self.propen_4[i + 1] = propensities[3]
            self.propen_5[i + 1] = propensities[4]
            self.propen_6[i + 1] = propensities[5]

    def update_particles(self, reaction, i):
        self.A[i + 1] = self.A[i]
        self.B[i + 1] = self.B[i]
        self.Aj[i + 1] = self.Aj[i]
        self.Bj[i + 1] = self.Bj[i]
        self.C[i + 1] = self.C[i]
        self.Cj[i + 1] = self.Cj[i]
        self.P[i + 1] = self.P[i]
        self.D[i + 1] = self.D[i]

        if reaction == 0:
            self.A[i + 1] -= 1
            self.B[i + 1] -= 1
            self.Aj[i + 1] += 1
        elif reaction == 1:
            self.Aj[i + 1] -= 1
            self.Bj[i + 1] += 1
            self.C[i + 1] -= 1
        elif reaction == 2:
            self.B[i + 1] += 1
            self.Bj[i + 1] -= 1
            self.C[i + 1] -= 2
            self.P[i + 1] += 1
        elif reaction == 3:
            self.A[i + 1] += 1
            self.B[i + 1] -= 1
            self.Cj[i + 1] += 1
            self.P[i + 1] -= 1
        elif reaction == 4:
            self.A[i + 1] -= 1
            self.B[i + 1] += 1
            self.Cj[i + 1] -= 1
            self.P[i + 1] += 1
        elif reaction == 5:
            self.Cj[i + 1] -= 1
            self.P[i + 1] -= 1
            self.D[i + 1] += 1

    def organize(self):
        rows = []
        for _ in range(self.params['Nc']):
            self.initialize_arrays()
            self.run_simulation()
            data = pd.DataFrame({
                'AParticles': [self.A], 'BParticles': [self.B], 'CParticles': [self.C],
                'PParticles': [self.P], 'AjParticles': [self.Aj], 'BjParticles': [self.Bj],
                'CjParticles': [self.Cj], 'DParticles': [self.D], 'Dif': [self.dif],
                'time': [self.time], 'Temperature': [self.T], 'sum_propen': [self.sum_propen],
                'propen_1': [self.propen_1], 'propen_2': [self.propen_2],
                'propen_3': [self.propen_3], 'propen_4': [self.propen_4],
                'propen_5': [self.propen_5], 'propen_6': [self.propen_6]
            })
            rows.append(data)
        return pd.concat(rows, ignore_index=True)
