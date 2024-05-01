import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as sci
from dataclasses import dataclass
from Planet import *

@dataclass
class Booster:
    thrust: float = 167970.0 #Newtons
    ISP1: float = 250.0 #s
    ISP2: float = 400.0 #s
    mass_stage1:float = 4300.0 #kg
    mass_stage2:float = 100.0 #kg

@dataclass
class Satellite:
    x: float
    vx: float
    z: float
    vz: float
    mass:float = 100.0 #kg
    
@dataclass
class Body:
    name: str
    booster: Booster
    satellite: Satellite
    altitud: float = 0.0

    def __post_init__(self):
        self.mass = self.satellite.mass + self.booster.mass_stage1 + self.booster.mass_stage2
        self.state = np.array([self.satellite.x, self.satellite.vx, self.satellite.z, self.satellite.vz, self.mass])
        self.r0 = np.sqrt(self.state[0]**2+self.state[1]**2)
        
@dataclass
class MotionSimulator:
    planet: Planet
    body: Body
    G: float = 6.6742e-11
    t_burn_1: float = 38.0
    t_step1: float = 40.0
    t_burn_2: float = 150.0
    t_end: float = 167.0
    def __post_init__(self):
        self.period = 2*np.pi/np.sqrt(self.G*self.planet.mass)*self.body.r0**(3/2)
        
    def gravity(self, x, z):
        r = np.sqrt(x**2 + z**2)
        ax = self.G*self.planet.mass/(r**3)*x
        az = self.G*self.planet.mass/(r**3)*z
        if r<0:
            ax, az = 0.0, 0.0
        return np.array([ax, az])
    
    def thrust(self, t):
        theta = 0.0
        if t<self.t_burn_1:
            theta = 10*np.pi/180
            thrust_F = self.body.booster.thrust
            self.ve = self.body.booster.ISP1*9.81
            mdot = -thrust_F/self.ve
        if (t>self.t_burn_1) and (t<self.t_step1):
            theta = 0.0
            thrust_F = 0.0
            mdot = -self.body.booster.mass_stage1/self.t_step1
        if t>self.t_step1:
            theta = 0.0
            thrust_F = 0.0
            mdot = 0.0
        if (t>self.t_burn_2) and (t<self.t_end):    
            theta = 90*np.pi/180
            thrust_F = self.body.booster.thrust
            self.ve = self.body.booster.ISP2*9.81
            mdot = -thrust_F/self.ve
        if t>self.t_end:
            theta = 0.0
            thrust_F = 0.0
            mdot = 0.0
            
            
        thrust_x = thrust_F*np.cos(theta)
        thrust_y = thrust_F*np.sin(theta)
        return np.array([thrust_x, thrust_y]), mdot
    
    def derivatives(self, state_in, t):
        x, xdot, z, zdot, mass = state_in
        
        gravity_F = -self.gravity(x, z) * mass
        thrust_F, mdot = self.thrust(t)
        aero_F = np.array([0.0, 0.0])
        
        F = gravity_F + aero_F + thrust_F
        
        if mass>0:
            ddot = F / mass
        else:
            ddot = [0.0, 0.0]
            mdot = 0.0

        state_out = [xdot, ddot[0], zdot, ddot[1], mdot]
        
        return state_out

    def simulate(self, t_total):
        state_out = sci.odeint(self.derivatives, self.body.state, t_total)
        return state_out

# Example usage:
planet = Planet("Kerbin")
booster = Booster()
satellite = Satellite(planet.radius, 0.0, 0.0, 0.0)
shuttle = Body("Vroom", booster, satellite) 
simulator = MotionSimulator(planet, shuttle)
t_total = np.linspace(0, 5000, 1000)
theta = np.linspace(0,2*np.pi,1000)
state_out = simulator.simulate(t_total)

plt.figure()
plt.plot(state_out[:, 0], state_out[:, 2])
plt.plot(planet.radius*np.sin(theta), planet.radius*np.cos(theta))
plt.xlabel('Time')
plt.ylabel('Position')
plt.title('Position vs Time')
plt.tight_layout()
plt.grid()
plt.show()

plt.figure()
plt.plot(t_total, state_out[:, 4])
plt.tight_layout()
plt.grid()
plt.show()

plt.figure()
plt.plot(t_total, state_out[:, 0])
plt.tight_layout()
plt.grid()
plt.show()