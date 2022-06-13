from FESTIM import Stepsize


class CyclingStepsize(Stepsize):
    def __init__(
        self, rampup, plateau, rampdown, rest, nb_cycles, stepsizes_max, **kwargs
    ) -> None:
        super().__init__(**kwargs)
        self.rampup = rampup
        self.rampdown = rampdown
        self.plateau = plateau
        self.rest = rest
        self.nb_cycles = nb_cycles
        self.stepsizes_max = stepsizes_max

    def adapt(self, t, nb_it, converged):
        """Changes the stepsize based on convergence.

        Args:
            t (float): current time.
            nb_it (int): number of iterations the solver required to converge.
            converged (bool): True if the solver converged, else False.
        """
        change_ratio = self.adaptive_stepsize["stepsize_change_ratio"]
        dt_min = self.adaptive_stepsize["dt_min"]
        stepsize_stop_max = self.adaptive_stepsize["stepsize_stop_max"]
        # t_stop = self.adaptive_stepsize["t_stop"]
        if not converged:
            self.value.assign(float(self.value) / change_ratio)
            if float(self.value) < dt_min:
                raise ValueError("stepsize reached minimal value")
        if nb_it < 5:
            self.value.assign(float(self.value) * change_ratio)
        else:
            self.value.assign(float(self.value) / change_ratio)

        if self.during_rampup(t) and float(self.value) > self.stepsizes_max["rampup"]:
            self.value.assign(self.stepsizes_max["rampup"])
        if self.during_plateau(t) and float(self.value) > self.stepsizes_max["plateau"]:
            self.value.assign(self.stepsizes_max["plateau"])
        if (
            self.during_rampdown(t)
            and float(self.value) > self.stepsizes_max["rampdown"]
        ):
            self.value.assign(self.stepsizes_max["rampdown"])
        if self.during_rest(t) and float(self.value) > self.stepsizes_max["rest"]:
            self.value.assign(self.stepsizes_max["rest"])

    def during_rampup(self, t):
        length_cycle = self.rampup + self.plateau + self.rampdown + self.rest
        current_cycle = t // length_cycle
        beginning_cycle = current_cycle * length_cycle
        return beginning_cycle <= t < beginning_cycle + self.rampup

    def during_plateau(self, t):
        length_cycle = self.rampup + self.plateau + self.rampdown + self.rest
        current_cycle = t // length_cycle
        beginning_cycle = current_cycle * length_cycle
        beggining_plateau = beginning_cycle + self.rampup
        return beggining_plateau <= t < beggining_plateau + self.plateau

    def during_rampdown(self, t):
        length_cycle = self.rampup + self.plateau + self.rampdown + self.rest
        current_cycle = t // length_cycle
        beginning_cycle = current_cycle * length_cycle
        beggining_rampdown = beginning_cycle + self.rampup + self.plateau
        return beggining_rampdown <= t < beggining_rampdown + self.rampdown

    def during_rest(self, t):
        return not (
            self.during_rampup(t) or self.during_plateau(t) or self.during_rampdown(t)
        )


if __name__ == "__main__":
    my_stepsize = CyclingStepsize(2, 10, 3, 50, nb_cycles=2, stepsizes_max=0)
    print(my_stepsize.during_rampup(66))
    print(my_stepsize.during_plateau(70))
    print(my_stepsize.during_rampdown(12))
    print(my_stepsize.during_rest(15))
