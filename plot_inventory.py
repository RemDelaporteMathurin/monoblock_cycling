import numpy as np
import matplotlib.pyplot as plt


data_cycling = np.genfromtxt(
    "results/cycling/derived_quantities.csv", delimiter=",", names=True
)


inventory_cycling = sum(
    [data_cycling["Total_retention_volume_{}".format(i)] for i in [1, 2, 3]]
)

data_continuous = np.genfromtxt(
    "results/continuous/derived_quantities.csv", delimiter=",", names=True
)
inventory_continuous = sum(
    [data_continuous["Total_retention_volume_{}".format(i)] for i in [1, 2, 3]]
)

plt.plot(data_cycling["ts"], inventory_cycling)
plt.plot(data_continuous["ts"], inventory_continuous)
# plt.xscale("log")
plt.yscale("log")
plt.show()
