from main import model, tungsten, recombination_flux_coolant, convection_flux

from cycling_stepsize import CyclingStepsize
from fluxes import CyclingFlux, CyclingImplantationDirichlet

import FESTIM as F

nb_cycles = 20
rampup = 100
plateau = 400
rampdown = 100
rest = 1000

heat_flux_value = 5e6
part_flux_value = 5e21

model.settings = F.Settings(
    absolute_tolerance=1e10,
    relative_tolerance=1e-10,
    final_time=nb_cycles * (rampup + plateau + rampdown + rest),
    chemical_pot=False,
    traps_element_type="DG",
)

model.dt = CyclingStepsize(
    rampup=rampup,
    plateau=plateau,
    rampdown=rampdown,
    rest=rest,
    nb_cycles=nb_cycles,
    stepsizes_max={"rampup": 5, "plateau": 20, "rampdown": 5, "rest": 50},
    initial_value=0.1,
    stepsize_change_ratio=1.1,
    dt_min=0.1,
)

h_implantation = CyclingImplantationDirichlet(
    surfaces=1,
    phi=part_flux_value,
    R_p=9.52e-10,
    D_0=tungsten.D_0,
    E_D=tungsten.E_D,
    data_y=[0, part_flux_value, part_flux_value, 0],
    rampup=rampup,
    plateau=plateau,
    rampdown=rampdown,
    rest=rest,
    nb_cycles=nb_cycles,
)

heat_flux = CyclingFlux(
    surfaces=1,
    field="T",
    data_y=[0, heat_flux_value, heat_flux_value, 0],
    rampup=rampup,
    plateau=plateau,
    rampdown=rampdown,
    rest=rest,
    nb_cycles=nb_cycles,
)


h_transport_bcs = [h_implantation, recombination_flux_coolant]
heat_transfer_bcs = [heat_flux, convection_flux]


model.boundary_conditions = h_transport_bcs + heat_transfer_bcs

model.T = F.HeatTransferProblem(
    transient=True, initial_value=323, relative_tolerance=1e-6, absolute_tolerance=1e0
)
main_folder = "results/phi_heat={:.1e}_phi_part={:.1e}/cycling".format(
    heat_flux_value, part_flux_value
)
derived_quantities = F.DerivedQuantities(
    [
        F.TotalVolume(field="retention", volume=1),
        F.TotalVolume(field="retention", volume=2),
        F.TotalVolume(field="retention", volume=3),
    ],
    filename=main_folder + "/derived_quantities.csv",
    nb_iterations_between_exports=10,
)

model.exports = F.Exports(
    [
        F.XDMFExport(
            "solute",
            filename=main_folder + "/mobile_concentration.xdmf",
            checkpoint=False,
        ),
        F.XDMFExport(
            "retention",
            filename=main_folder + "/retention.xdmf",
            checkpoint=False,
        ),
        F.XDMFExport("T", filename=main_folder + "/temperature.xdmf", checkpoint=False),
        derived_quantities,
    ]
)
model.log_level = 20
model.initialise()
model.run()
