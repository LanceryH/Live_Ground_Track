import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.basemap import Basemap
import scipy.integrate as sci
from dataclasses import dataclass
from datetime import datetime, timezone
import plotly.graph_objs as go
from Planet import *
from Display import *
from Body import *
from Simulator import *
def conversion_to_geodic(x,y,z,centered=False):
        #if len(R.shape)==1:
        #    R=R.reshape((3,1))
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
  
# Example usage:
planet = Planet("Kerbin")
booster = Booster()
satellite = Satellite(planet.radius, 0.0, 0.0, 0.0, 0.0, 0.0)
shuttle = Body("Vroom", booster, satellite) 
simulator = MotionSimulator(planet, shuttle)
t_total = np.linspace(0, 5000, 1000)
theta = np.linspace(0,2*np.pi,1000)
state_out = simulator.simulate(t_total)


lon, lat = conversion_to_geodic(state_out[:, 0],state_out[:, 2],state_out[:, 4], centered=True)

#np.savetxt("satellite.csv", [lat,lon], delimiter=",")

#lat[np.abs(lon)+180 <= lon[1]-lon[0]] = np.nan
lon_plot = lon.tolist()
lat_plot = lat.tolist()

for ind_i in range(0,len(lon_plot)-1):
    if np.abs(lon_plot[ind_i]+lon_plot[ind_i+1])<np.abs(lon_plot[ind_i]):
        lon_plot[ind_i]=np.nan
        lat_plot[ind_i]=np.nan

# miller projection 
plt.figure()
map = Basemap(projection='mill',lon_0=0)
map.drawcoastlines(color="yellow")
map.drawparallels(np.arange(-90,90,30),labels=[1,0,0,0])
map.drawmeridians(np.arange(map.lonmin,map.lonmax+30,60),labels=[0,0,0,1])
# fill continents 'coral' (with zorder=0), color wet areas 'aqua'
map.drawmapboundary(fill_color='lightsteelblue')
map.fillcontinents(color='darkseagreen',lake_color='lightsteelblue')
map.drawcountries(color="yellow")
# shade the night areas, with alpha transparency so the 
# map shows through. Use current time in UTC.
CS=map.nightshade(datetime.now(timezone.utc))
x, y = map(lon_plot, lat_plot)
x_0, y_0 = map(lon_plot[0], lat_plot[0])
plt.plot(x, y, color="r")
plt.plot(x_0, y_0, color="r", marker="o")
plt.show()