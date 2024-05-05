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
    t_burn_1: float = 146
    t_step1: float = 2.0
    t_burn_2: float = 148
    t_end: float = 148+138
    
    def __post_init__(self):
        self.period = 2*np.pi/np.sqrt(self.G*self.planet.mass)*self.body.r0**(3/2)
    
        
    def gravity_func(self, x, y, z, r):
        ax = self.G*self.planet.mass/(r**3)*x
        ay = self.G*self.planet.mass/(r**3)*y
        az = self.G*self.planet.mass/(r**3)*z
        if r<self.planet.radius:
            ax, ay, az = 0.0, 0.0, 0.0
        return np.array([ax, ay, az])
    
    
    def temp_rho_func(self, altitud):
        if self.planet.name == "Kerbin":
            temp = 25
            pressure = 1000
            rho = np.interp(altitud,self.planet.altitud, self.planet.density)
        if self.planet.name == "Earth":
            if altitud < 11000:
                temp = 15.04 - 0.00649*altitud
                pressure = 101.29 * ((temp + 273.15)/288.08)**5.256
            if (altitud >11000) & (altitud<25000):
                temp = -56.49
                pressure = 22.65*np.exp(1.73-0.000157*altitud)
            if altitud>25000:
                temp = -131.21+0.00299* altitud
                pressure = 2.488 * ((temp + 273.15)/216.6)**-11.388
            rho = pressure/(0.2869*(temp+273.15))
            
        return temp, pressure, rho
    
    def thrust_func(self, t, m, x, y, z, r):
        theta_live = np.arccos(z/r)
        phi_live = np.arctan(y/x)
        theta = theta_live
        phi = phi_live
        
        if t<self.t_burn_1:
            #theta = np.deg2rad(80)
            #phi = np.deg2rad(45)
            theta = np.deg2rad(80)
            phi = np.deg2rad(45)
            if m>(self.body.mass-self.body.booster.mass_stage1):
                thrust_F = self.body.booster.thrust_stage_1
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
            #theta = np.deg2rad(70)#+theta_live
            #phi = np.deg2rad(120)#-phi_live
            theta = np.deg2rad(80)
            phi = np.deg2rad(45)
            if m>(self.body.mass-self.body.booster.mass_stage1-self.body.booster.mass_stage2):
                thrust_F = self.body.booster.thrust_stage_2
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
    
    def aero_func(self, v, rho, vx, vy, vz):
        return -0.5*self.body.surface*rho*self.body.drag_coef*v*np.array([vx, vy, vz])
    
    def derivatives(self, state_in, t):
        x, xdot, y, ydot, z, zdot, r_old, v_old, mass, temperatur_old, pressure_old, rho_old  = state_in
        r = np.sqrt(x**2+y**2+z**2)
        v = np.sqrt(xdot**2+ydot**2+zdot**2)
        altitud = r-self.planet.radius
        temperatur, pressure, rho = self.temp_rho_func(altitud)
        gravity_F = -self.gravity_func(x, y, z, r) * mass
        thrust_F, mdot = self.thrust_func(t, mass, x, y, z, r)
        aero_F = self.aero_func(v, rho, xdot, ydot, zdot)
        F = gravity_F + aero_F + thrust_F
            
        if r<self.planet.radius:
            return np.zeros_like(state_in)
        else:
            ddot = F / mass
            return [xdot, ddot[0], ydot, ddot[1], zdot, ddot[2], r-r_old, v-v_old, mdot, 
                    temperatur-temperatur_old, 
                    pressure-pressure_old, 
                    rho-rho_old]        

    def simulate(self, t_total):
        self.state_out = sci.odeint(self.derivatives, self.body.state, t_total)
        return self.state_out
    
    def cart2geodic(self,centered=False):
        x, y, z = self.state_out[:,0::2]
        p=np.sqrt(x**2+y**2+z**2)
        #GMST = self.calculate_GMST()
        lon=np.rad2deg(np.arctan2(y,x))#-GMST
        lat=np.rad2deg(np.arcsin(z/p))
        for j in range(len(lon)):
            if lon[j] > 180:
                lon[j] -= 360
            elif lon[j] < -180:
                lon[j] += 360
            #a = 6378137
            #b = 6356752.314235
            #f = (a - b) / a
            #lat[j] = lat[j] * (1 - f * f)
        return lon, lat