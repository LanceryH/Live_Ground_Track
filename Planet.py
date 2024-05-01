from dataclasses import dataclass

@dataclass
class Planet:
    name: str
    mass: float = 0.0
    radius: float = 0.0
    
    def __post_init__(self):
        if self.name == "Earth":
            self.mass = 5.9722e24
            self.radius = 6378000
            self.gravity = 9.81
        if self.name == "Kerbin":
            self.mass = 5.2982e22
            self.radius = 6e5
            self.gravity = 9.81