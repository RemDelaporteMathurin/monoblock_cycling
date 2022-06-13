import fenics as f
from scipy.interpolate import interp1d
import sympy as sp
import FESTIM as F


class InterpolatedExpression(f.UserExpression):
    def __init__(self, data_t, data_y) -> None:
        super().__init__()
        self.interpolated_object = interp1d(data_t, data_y, fill_value=0)

    def eval(self, value, x):
        value[0] = self.interpolated_object(self.t)


class CyclingFlux(F.FluxBC):
    def __init__(
        self, rampup, plateau, rampdown, rest, nb_cycles, data_y, **kwargs
    ) -> None:
        super().__init__(**kwargs)
        self.rampup = rampup
        self.rampdown = rampdown
        self.plateau = plateau
        self.rest = rest
        self.nb_cycles = nb_cycles
        assert len(data_y) == 4
        self.data_y = data_y * nb_cycles + [data_y[-1]]

    def create_form(self, T, solute):
        self.form = InterpolatedExpression(
            data_t=self.make_t_data(), data_y=self.data_y
        )
        self.sub_expressions.append(self.form)

    def make_t_data(self):
        data_t = []
        cycle_length = self.rampup + self.plateau + self.rampdown + self.rest
        for n in range(self.nb_cycles):
            data_t += [
                n * cycle_length,
                self.rampup + n * cycle_length,
                (self.rampup + self.plateau) + n * cycle_length,
                (self.rampup + self.plateau + self.rampdown) + n * cycle_length,
            ]

        data_t.append(data_t[-1] + self.rest)

        return data_t


class CyclingImplantationDirichlet(F.ImplantationDirichlet):
    def __init__(
        self, rampup, plateau, rampdown, rest, nb_cycles, data_y, **kwargs
    ) -> None:
        super().__init__(**kwargs)
        self.rampup = rampup
        self.rampdown = rampdown
        self.plateau = plateau
        self.rest = rest
        self.nb_cycles = nb_cycles
        assert len(data_y) == 4
        self.data_y = data_y * nb_cycles + [data_y[-1]]

    def create_expression(self, T):
        phi = InterpolatedExpression(data_t=self.make_t_data(), data_y=self.data_y)
        R_p = f.Expression(sp.printing.ccode(self.R_p), t=0, degree=1)
        sub_expressions = [phi, R_p]

        value_BC = F.BoundaryConditionExpression(
            T,
            dc_imp,
            phi=phi,
            R_p=R_p,
            D_0=self.D_0,
            E_D=self.E_D,
            Kr_0=self.Kr_0,
            E_Kr=self.E_Kr,
        )
        self.expression = value_BC
        self.sub_expressions = sub_expressions

    def make_t_data(self):
        data_t = []
        cycle_length = self.rampup + self.plateau + self.rampdown + self.rest
        for n in range(self.nb_cycles):
            data_t += [
                n * cycle_length,
                self.rampup + n * cycle_length,
                (self.rampup + self.plateau) + n * cycle_length,
                (self.rampup + self.plateau + self.rampdown) + n * cycle_length,
            ]
        data_t.append(data_t[-1] + self.rest)

        return data_t


def dc_imp(T, phi, R_p, D_0, E_D, Kr_0=None, E_Kr=None):
    D = D_0 * f.exp(-E_D / F.k_B / T)
    value = phi * R_p / D
    if Kr_0 is not None:
        Kr = Kr_0 * f.exp(-E_Kr / F.k_B / T)
        value += (phi / Kr) ** 0.5

    return value


# data_particle_flux = [0, 1.61e22, 1.61e22, 0]

# heat_flux = InterpolatedExpression(data_t=data_t, data_y=data_heat_flux)
# particle_flux = InterpolatedExpression(data_t=data_t, data_y=data_particle_flux)
