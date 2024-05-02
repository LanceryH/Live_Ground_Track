from dataclasses import dataclass

@dataclass
class Booster:
    thrust: float = 150000.0 #Newtons
    ISP1: float = 250.0 #s
    ISP2: float = 400.0 #s
    mass_stage1:float = 4000.0 #kg
    mass_stage2:float = 1000.0 #kg