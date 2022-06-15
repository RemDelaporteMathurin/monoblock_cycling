import numpy as np
import matplotlib.pyplot as plt
import matplotx
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit

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


def get_inventory_cycling(heat_flux, part_flux):
    folder = "results/phi_heat={:.1e}_phi_part={:.1e}".format(heat_flux, part_flux)
    data_cycling = np.genfromtxt(
        folder + "/cycling/derived_quantities.csv", delimiter=",", names=True
    )

    inventory_cycling = sum(
        [data_cycling["Total_retention_volume_{}".format(i)] for i in [1, 2, 3]]
    )
    flux = get_flux(part_flux)

    dt_cycling = np.diff(data_cycling["ts"])
    fluence_cycling = [0]
    for i, t in enumerate(data_cycling["ts"][:-1]):
        fluence_cycling.append(fluence_cycling[-1] + flux(t) * dt_cycling[i])

    fluence_cycling = np.array(fluence_cycling)
    return fluence_cycling, inventory_cycling


def get_inventory_continuous(heat_flux, part_flux):
    folder = "results/phi_heat={:.1e}_phi_part={:.1e}".format(heat_flux, part_flux)
    data_continuous = np.genfromtxt(
        folder + "/continuous/derived_quantities.csv", delimiter=",", names=True
    )
    inventory_continuous = sum(
        [data_continuous["Total_retention_volume_{}".format(i)] for i in [1, 2, 3]]
    )

    fluence_continuous = data_continuous["ts"] * part_flux
    return fluence_continuous, inventory_continuous


def plot_continuous_cycling(heat_flux, part_flux, label="", **kwargs):

    fluence_cycling, inventory_cycling = get_inventory_cycling(heat_flux, part_flux)
    fluence_continuous, inventory_continuous = get_inventory_continuous(
        heat_flux, part_flux
    )

    plt.plot(
        fluence_cycling[np.where(fluence_cycling < fluence_max)],
        inventory_cycling[np.where(fluence_cycling < fluence_max)],
        **kwargs
    )

    plt.plot(
        fluence_continuous[np.where(fluence_continuous < fluence_max)],
        inventory_continuous[np.where(fluence_continuous < fluence_max)],
        linestyle="dashed",
        label=label,
        **kwargs
    )


def fit_continuous(heat_flux, part_flux):
    def power_law(x, a, b):
        return a * x**b

    folder = "results/phi_heat={:.1e}_phi_part={:.1e}".format(heat_flux, part_flux)
    data_continuous = np.genfromtxt(
        folder + "/continuous/derived_quantities.csv", delimiter=",", names=True
    )
    inventory_continuous = sum(
        [data_continuous["Total_retention_volume_{}".format(i)] for i in [1, 2, 3]]
    )
    fluence_continuous = data_continuous["ts"] * part_flux

    popt, pcov = curve_fit(power_law, fluence_continuous, inventory_continuous)
    return popt


fluence_max = 3.6e25


plt.figure(figsize=(6.4, 3))

plot_continuous_cycling(5e6, 5e21, color="tab:blue", label="low flux")
a, b = fit_continuous(5e6, 5e21)
plt.annotate(
    "$\propto \mathrm{fluence} ^{" + "{:.1f}".format(b) + "}$",
    (1e24, 5e19),
    color="tab:blue",
)


plot_continuous_cycling(13e6, 1.6e22, color="tab:orange", label="high flux")
a, b = fit_continuous(13e6, 1.6e22)
plt.annotate(
    "$\propto \mathrm{fluence} ^{" + "{:.1f}".format(b) + "}$",
    (2e24, 3e18),
    color="tab:orange",
)

plt.xlim(left=6e23)
plt.ylim(4e17, 7e20)
plt.xscale("log")
plt.yscale("log")

plt.gca().spines.right.set_visible(False)
plt.gca().spines.top.set_visible(False)
plt.xlabel("Fluence (m$^{-2}$)")
plt.ylabel("Inventory (m$^{-2}$)")
matplotx.line_labels()
plt.tight_layout()
plt.show()
