from dataclasses import dataclass
import numpy as np

@dataclass
class Planet:
    name: str
    
    def __post_init__(self):
        if self.name == "Earth":
            self.mass = 5.9722e24
            self.radius = 6378000
            self.gravity = 9.81
            self.beta = 0.1354/1000.0 # density cte
            self.rhos = 1.225 # kg/m³
        if self.name == "Kerbin":
            self.mass = 5.2982e22
            self.radius = 6e5
            self.gravity = 9.81
            self.atmo = np.loadtxt("kerbin_atmo.txt")
            self.altitud = self.atmo[:,0]
            self.density = self.atmo[:,3]
            self.rhos = self.density[0]
            self.beta = 0.0
            
            
