import FESTIM as F

import numpy as np


model = F.Simulation()


def thermal_cond_W(T):
    return -7.84154e-9 * T**3 + 5.03006e-5 * T**2 - 1.07335e-1 * T + 1.75214e2


def thermal_cond_Cu(T):
    return -3.93153e-08 * T**3 + 3.76147e-05 * T**2 - 7.88669e-02 * T + 4.02301e02


def thermal_cond_CuCrZr(T):
    return 5.25780e-7 * T**3 - 6.45110e-4 * T**2 + 2.57678e-01 * T + 3.12969e2


def rho_cp_W(T):
    return 5.15356e-6 * T**3 - 8.30703e-2 * T**2 + 5.98312e2 * T + 2.48160e6


def rho_cp_Cu(T):
    return 1.68402e-4 * T**3 - 6.14079e-2 * T**2 + 4.67353e2 * T + 3.45899e6


def rho_cp_CuCrZr(T):
    return -1.79134e-4 * T**3 - 1.51383e-1 * T**2 + 6.22091e2 * T + 3.46007e6


# atom_density  =  density(g/m3)*Na(/mol)/M(g/mol)
atom_density_W = 6.28e28  # 6.3222e28  # atomic density m^-3
atom_density_Cu = 8.43e28  # 8.4912e28  # atomic density m^-3
atom_density_CuCrZr = 2.6096e28  # atomic density m^-3

tungsten = F.Material(
    id=1,
    D_0=1.9e-7,
    E_D=0.2,
    S_0=1.87e24,
    E_S=1.04,
    thermal_cond=thermal_cond_W,
    rho=1,
    heat_capacity=rho_cp_W,
    borders=[0, 6e-3],
)
cu = F.Material(
    id=2,
    D_0=6.6e-7,
    E_D=0.39,
    S_0=3.14e24,
    E_S=0.57,
    thermal_cond=thermal_cond_Cu,
    rho=1,
    heat_capacity=rho_cp_Cu,
    borders=[6e-3, 7e-3],
)
cucrzr = F.Material(
    id=3,
    D_0=3.9e-7,
    E_D=0.42,
    S_0=4.28e23,
    E_S=0.39,
    thermal_cond=thermal_cond_CuCrZr,
    rho=1,
    heat_capacity=rho_cp_CuCrZr,
    borders=[7e-3, 7.5e-3],
)

model.materials = [tungsten, cu, cucrzr]

trap_w2 = F.Trap(
    8.96e-17, 0.2, 1e13, 1, materials=tungsten, density=4e-4 * atom_density_W
)
# trap_cu = F.Trap(
#     6.0e-17, 0.39, 8.0e13, 0.50, materials=cu, density=5.0e-5 * atom_density_Cu
# )
# trap_cucrzr = F.Trap(
#     1.2e-16, 0.42, 8.0e13, 0.85, materials=cucrzr, density=5.0e-5 * atom_density_CuCrZr
# )

trap_conglo = F.Trap(
    k_0=[8.96e-17, 6.0e-17, 1.2e-16],
    E_k=[0.2, 0.39, 0.42],
    p_0=[1e13, 8e13, 8e13],
    E_p=[0.87, 0.5, 0.85],
    materials=[tungsten, cu, cucrzr],
    density=[
        1.1e-3 * atom_density_W,
        5.0e-5 * atom_density_Cu,
        5.0e-5 * atom_density_CuCrZr,
    ],
)

model.traps = [trap_conglo, trap_w2]

# model.mesh = F.MeshFromVertices(
#     np.concatenate(
#         [
#             np.linspace(tungsten.borders[0], 1e-6, num=50),
#             np.linspace(1e-6, tungsten.borders[1], num=500),
#             np.linspace(*cu.borders, num=50),
#             np.linspace(*cucrzr.borders, num=100),
#         ]
#     )
# )

model.mesh = F.MeshFromRefinements(
    600,
    size=cucrzr.borders[-1],
    refinements=[
        {"x": tungsten.borders[-1], "cells": 500},
        {"x": 1e-5, "cells": 100},
        {"x": 2e-6, "cells": 200},
        {"x": 2e-7, "cells": 20},
        # {"x": 1e-4, "cells": 50},
        # {"x": 1e-4, "cells": 100},
        # {"x": 1e-6, "cells": 100},
        # {"x": 1e-7, "cells": 100},
    ],
)

convection_flux = F.ConvectiveFlux(h_coeff=70000, T_ext=323, surfaces=2)

recombination_flux_coolant = F.RecombinationFlux(
    Kr_0=2.9e-14, E_Kr=1.92, order=2, surfaces=2
)
