from dataclasses import dataclass

@dataclass
class Satellite:
    x: float
    vx: float
    y: float
    vy: float
    z: float
    vz: float
    mass:float = 1000.0 #kg
    