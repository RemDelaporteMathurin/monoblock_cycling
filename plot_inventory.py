import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

nb_cycles = 30
rampup = 100
plateau = 400
rampdown = 100
rest = 1000


def get_flux(part_flux_value):
    data_t = []
    cycle_length = rampup + plateau + rampdown + rest
    for n in range(nb_cycles):
        data_t += [
            n * cycle_length,
            rampup + n * cycle_length,
            (rampup + plateau) + n * cycle_length,
            (rampup + plateau + rampdown) + n * cycle_length,
        ]

    data_t.append(data_t[-1] + rest)

    data_flux_part = [0, part_flux_value, part_flux_value, 0]
    data_flux_part = data_flux_part * nb_cycles + [data_flux_part[-1]]

    flux = interp1d(data_t, data_flux_part)
    return flux


def plot_continuous_cycling(heat_flux, part_flux):

    folder = "results/phi_heat={:.1e}_phi_part={:.1e}".format(heat_flux, part_flux)
    data_cycling = np.genfromtxt(
        folder + "/cycling/derived_quantities.csv", delimiter=",", names=True
    )

    inventory_cycling = sum(
        [data_cycling["Total_retention_volume_{}".format(i)] for i in [1, 2, 3]]
    )

    data_continuous = np.genfromtxt(
        folder + "/continuous/derived_quantities.csv", delimiter=",", names=True
    )
    inventory_continuous = sum(
        [data_continuous["Total_retention_volume_{}".format(i)] for i in [1, 2, 3]]
    )

    flux = get_flux(part_flux)

    dt_cycling = np.diff(data_cycling["ts"])
    fluence_cycling = [0]
    for i, t in enumerate(data_cycling["ts"][:-1]):
        fluence_cycling.append(fluence_cycling[-1] + flux(t) * dt_cycling[i])

    plt.plot(fluence_cycling, inventory_cycling)

    fluence_continuous = data_continuous["ts"] * part_flux
    plt.plot(fluence_continuous, inventory_continuous)


for heat_flux, part_flux in zip([5e6, 13e6], [5e21, 1.6e22]):
    plot_continuous_cycling(heat_flux, part_flux)

plt.xlim(left=6e23, right=3.6e25)
plt.ylim(4e17, 7e20)
plt.xscale("log")
plt.yscale("log")

plt.xlabel("Fluence (m$^{-2}$)")
plt.ylabel("Inventory (m$^{-2}$)")
plt.show()
