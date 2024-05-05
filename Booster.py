from dataclasses import dataclass

@dataclass
class Booster:
    thrust_stage_1: float = 245_000.0 #Newtons
    thrust_stage_2: float = 73_000.0 #Newtons
    ISP1: float = 262.0 #s
    ISP2: float = 262.0 #s
    mass_stage1:float = 13_270.0 #kg
    mass_stage2:float = 3_285.0 #kg