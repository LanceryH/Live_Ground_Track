from dataclasses import dataclass

@dataclass
class Satellite:
    x: float # m
    vx: float # m.s⁻¹
    y: float # m
    vy: float # m.s⁻¹ 
    z: float # m
    vz: float # m.s⁻¹ 
    mass: float = 100.0 # kg
    temperatur: float = 25 # °C
    rho: float = 1.125 # kg/m³
    pressure: float = 1013.25 # hPa
    