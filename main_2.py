import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as sci
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

@dataclass
class Body:
    name: str
    x: float
    vx: float
    z: float
    vz: float
    mass: float
    
    altitud: float = 0.0
    
    def __post_init__(self):
        self.state = np.array([self.x, self.vx, self.z, self.vz, self.mass])
        self.r0 = np.sqrt(self.state[0]**2+self.state[1]**2)

@dataclass
class Booster:
    thrust: float #Newtons
    ISP: float = 200 #s
        
@dataclass
class MotionSimulator:
    planet: Planet
    booster: Booster
    body: Body
    G: float = 6.6742e-11
    
    def __post_init__(self):
        self.period = 2*np.pi/np.sqrt(self.G*self.planet.mass)*self.body.r0**(3/2)
        
    def gravity(self, x, z):
        r = np.sqrt(x**2 + z**2)
        ax = self.G*self.planet.mass/(r**3)*x
        az = self.G*self.planet.mass/(r**3)*z
        return np.array([ax, az])
    
    def thrust(self, t):
        theta = 0.0
        if t<5:
            thrust_F = self.booster.thrust
            thrust_x = thrust_F*np.cos(theta)
            thrust_y = thrust_F*np.sin(theta)
        else:
            thrust_F = 0.0
            thrust_x = 0.0
            thrust_y = 0.0
        ve = self.booster.ISP*9.81
        mdot = -thrust_F/ve
        return np.array([thrust_x, thrust_y]), mdot
    
    def derivatives(self, state_in, t):
        x, xdot, z, zdot, mass = state_in
        
        gravity_F = -self.gravity(x, z) * mass
        thrust_F, mdot = self.thrust(t)
        aero_F = np.array([0.0, 0.0])
        
        F = gravity_F + aero_F + thrust_F
        ddot = F / mass
        
        state_out = [xdot, ddot[0], zdot, ddot[1], mdot]
        
        return state_out

    def simulate(self, t_total):
        state_out = sci.odeint(self.derivatives, self.body.state, t_total)
        return state_out

# Example usage:
earth = Planet("Earth")
booster = Booster(20)
v_orbt = np.sqrt(MotionSimulator.G*earth.mass/(earth.radius+6e5))
satellite = Body("Sat", earth.radius, 0.0, 0.0, 0.0, 0.64) #[x,vx,z,vz]
simulator = MotionSimulator(earth, booster, satellite)
t_total = np.linspace(0, 25, 1000)
theta = np.linspace(0,2*np.pi,1000)
state_out = simulator.simulate(t_total)

plt.figure()
plt.plot(state_out[:, 0], state_out[:, 2])
plt.plot(earth.radius*np.sin(theta), earth.radius*np.cos(theta))
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