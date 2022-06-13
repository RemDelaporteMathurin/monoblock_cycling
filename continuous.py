from main import model, tungsten, recombination_flux_coolant, convection_flux

import FESTIM as F

nb_cycles = 30
rampup = 100
plateau = 400
rampdown = 100
rest = 1000

model.settings = F.Settings(
    absolute_tolerance=1e10,
    relative_tolerance=1e-10,
    final_time=nb_cycles * (rampup + plateau + rampdown),
    chemical_pot=True,
    traps_element_type="DG",
)

model.dt = F.Stepsize(initial_value=1, stepsize_change_ratio=1.1, dt_min=0.1)

h_implantation = F.ImplantationDirichlet(
    surfaces=1, phi=1.61e22, R_p=9.52e-10, D_0=tungsten.D_0, E_D=tungsten.E_D
)

heat_flux = F.FluxBC(surfaces=1, field="T", value=10e6)


h_transport_bcs = [h_implantation, recombination_flux_coolant]
heat_transfer_bcs = [heat_flux, convection_flux]


model.boundary_conditions = h_transport_bcs + heat_transfer_bcs

model.T = F.HeatTransferProblem(transient=False)

derived_quantities = F.DerivedQuantities(
    [
        F.TotalVolume(field="retention", volume=1),
        F.TotalVolume(field="retention", volume=2),
        F.TotalVolume(field="retention", volume=3),
    ],
    filename="results/continuous/derived_quantities.csv",
)

model.exports = F.Exports(
    [
        F.XDMFExport(
            "solute",
            filename="results/continuous/mobile_concentration.xdmf",
            checkpoint=False,
        ),
        F.XDMFExport(
            "retention",
            filename="results/continuous/retention.xdmf",
            checkpoint=False,
        ),
        F.XDMFExport(
            "T", filename="results/continuous/temperature.xdmf", checkpoint=False
        ),
        derived_quantities,
    ]
)

model.initialise()
model.run()
