import numpy as np
from Planet import *
from Body import *
from dataclasses import dataclass
import scipy.integrate as sci

@dataclass
class MotionSimulator:
    planet: Planet
    body: Body
    G: float = 6.6742e-11
    t_burn_1: float = 63.0
    t_step1: float = 2.0
    t_burn_2: float = 250.0
    t_end: float = 263.0
    
    def __post_init__(self):
        self.period = 2*np.pi/np.sqrt(self.G*self.planet.mass)*self.body.r0**(3/2)
    
    def rho(self, x):
        if self.planet.name != "Eearth":
            rho = np.interp(x,self.planet.altitud, self.planet.density)
        else:
            rho = self.planet.rhos*np.exp(-self.planet.beta*x)
        return rho
        
    def gravity(self, x, y, z, r):
        ax = self.G*self.planet.mass/(r**3)*x
        ay = self.G*self.planet.mass/(r**3)*y
        az = self.G*self.planet.mass/(r**3)*z
        if r<self.planet.radius:
            ax, ay, az = 0.0, 0.0, 0.0
        return np.array([ax, ay, az])
    
    def thrust(self, t, m, x, y, z, r):
        theta_live = np.arccos(z/r)
        phi_live = np.arctan(y/x)
        theta = theta_live
        phi = phi_live
        
        if t<self.t_burn_1:
            theta = np.deg2rad(70)
            phi = np.deg2rad(45)
            if m>(self.body.mass-self.body.booster.mass_stage1):
                thrust_F = self.body.booster.thrust
                self.ve = self.body.booster.ISP1*9.81
                mdot = -thrust_F/self.ve
            else:
                thrust_F = 0.0
                mdot = 0.0
                
        if (t>self.t_burn_1) and (t<self.t_burn_1+self.t_step1):
            thrust_F = 0.0
            if m>(self.body.mass-self.body.booster.mass_stage1):
                mdot = -self.body.booster.mass_stage1/self.t_step1
            else:
                mdot = 0.0
                
        if t>(self.t_burn_1+self.t_step1):
            thrust_F = 0.0
            mdot = 0.0
            
        if (t>self.t_burn_2) and (t<self.t_end):    
            theta = np.deg2rad(70)#+theta_live
            phi = np.deg2rad(120)#-phi_live
            if m>(self.body.mass-self.body.booster.mass_stage1-self.body.booster.mass_stage2):
                thrust_F = self.body.booster.thrust
                self.ve = self.body.booster.ISP2*9.81
                mdot = -thrust_F/self.ve
            else:
                thrust_F = 0.0
                mdot = 0.0
                
        if t>self.t_end:
            thrust_F = 0.0
            mdot = 0.0
            
        thrust_x = thrust_F*np.sin(theta)*np.cos(phi)
        thrust_y = thrust_F*np.sin(theta)*np.sin(phi)
        thrust_z = thrust_F*np.cos(theta)
        return np.array([thrust_x, thrust_y, thrust_z]), mdot
    
    def aero(self, v, rho, vx, vy, vz):
        return -0.5*self.body.surface*rho*self.body.drag_coef*v*np.array([vx, vy, vz])
    
    def derivatives(self, state_in, t):
        x, xdot, y, ydot, z, zdot, mass = state_in
        r = np.sqrt(x**2+y**2+z**2)
        v = np.sqrt(xdot**2+ydot**2+zdot**2)
        altitud = r-self.planet.radius
        rho = self.rho(altitud)
        gravity_F = -self.gravity(x, y, z, r) * mass
        thrust_F, mdot = self.thrust(t, mass, x, y, z, r)
        aero_F = self.aero(v, rho, xdot, ydot, zdot)
        
        F = gravity_F + aero_F + thrust_F
            
        if r<self.planet.radius:
            return [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            
        else:
            ddot = F / mass
            return [xdot, ddot[0], ydot, ddot[1], zdot, ddot[2], mdot]        

    def simulate(self, t_total):
        state_out = sci.odeint(self.derivatives, self.body.state, t_total)
        return state_out