
class Fluid:
    def __init__(self,
                 density: float,
                 viscosity: float,
                 thermal_conductivity: float,
                 specific_heat_capacity: float):
        self.rho: float = density
        self.mu: float = viscosity
        self.alpha: float = thermal_conductivity
        self.cp: float = specific_heat_capacity

    def __str__(self) -> str:
        return (
            f'Density ρ, [kg/m^3] {self.rho}' +
            f'Dynamic viscosity μ,[Pa*s] {self.mu} ' +
            f'Thermal conductivity α, [W/(m*K)] {self.alpha}' +
            f'Specific heat capacity Cp, [J/kg*K] {self.cp}')
