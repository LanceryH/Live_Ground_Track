from Booster import *
from Satellite import *
from dataclasses import dataclass
import numpy as np

@dataclass
class Body:
    name: str
    booster: Booster
    satellite: Satellite
    radius: float = 0.5 #m
    height: float = 1 #m
    drag_coef: float = 0.4225
    
    def __post_init__(self):
        self.surface = np.pi*self.radius**2*self.height
        self.mass = self.satellite.mass + self.booster.mass_stage1 + self.booster.mass_stage2
        self.state = np.array([self.satellite.x, self.satellite.vx, self.satellite.y, self.satellite.vy, self.satellite.z, self.satellite.vz, self.mass])
        self.r0 = np.sqrt(self.state[0]**2+self.state[1]**2)
