import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime, timezone
import pytz
import pandas as pd
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import matplotlib.patches as ptc
import time

class Satellite_class:
    def __init__(self, data):
        self.MUE = 398600.44
        self.EARTH_MASS = 5.972e24
        self.G = 6.67384e-11
        self.NAME = data[0]
        self.NORA_ID = data[1]
        self.EPOCH = data[2]
        self.MEAN_MOTION = data[3]
        self.ECCENTRICITY = data[4]*np.pi/180
        self.INCLINATION = data[5]*np.pi/180
        self.RA_OF_ASC_NODE = data[6]*np.pi/180
        self.ARG_OF_PERICENTER = data[7]*np.pi/180
        self.MEAN_ANOMALY = data[8]*np.pi/180
        self.UPDATE_DATE_YEAR=int(self.EPOCH[:4])
        self.UPDATE_DATE_MONTH=int(self.EPOCH[5:7])
        self.UPDATE_DATE_DAY=int(self.EPOCH[8:10])
        self.UPDATE_DATE_HOUR=int(self.EPOCH[11:13])
        self.UPDATE_DATE_MINUTE=int(self.EPOCH[14:16])
        self.UPDATE_DATE_SECOND=int(self.EPOCH[17:19])
        self.MEAN_ANOMALY_UPDATED = None
        self.SEMI_MAJOR_AXIS = None
        self.ECCENTRIC_ANOMALY = None
        self.R = None
        self.Rdot = None
        self.lon = None
        self.lat = None

    def mean_motion(self):
        return np.sqrt(self.EARTH_MASS*self.G/self.SEMI_MAJOR_AXIS**3)

    def mean_anomaly(self,t):
        return self.MEAN_ANOMALY+self.MEAN_MOTION*t

    def Ecc(self):
        return self.MEAN_ANOMALY_UPDATED-self.ECCENTRICITY*np.sin(self.MEAN_ANOMALY_UPDATED)-self.MEAN_ANOMALY_UPDATED

    def Ecc_dot(self):
        return 1-self.ECCENTRICITY*np.cos(self.MEAN_ANOMALY_UPDATED)

    def eccentric_anomaly(self):
        epsilon = 1e-3
        convergence = 1
        ECCENTRIC_ANOMALY_N = self.MEAN_ANOMALY_UPDATED

        while convergence > epsilon:
            numerator = (ECCENTRIC_ANOMALY_N - self.ECCENTRICITY * np.sin(ECCENTRIC_ANOMALY_N) - self.MEAN_ANOMALY_UPDATED)
            denominator = (1 - self.ECCENTRICITY * np.cos(ECCENTRIC_ANOMALY_N))
            ECCENTRIC_ANOMALY_N1 = ECCENTRIC_ANOMALY_N + numerator / denominator
            convergence = np.abs(ECCENTRIC_ANOMALY_N1 - ECCENTRIC_ANOMALY_N)
            ECCENTRIC_ANOMALY_N = ECCENTRIC_ANOMALY_N1

        return ECCENTRIC_ANOMALY_N

    def R_INCLINATION(self):
        return np.array([[1,0,0],
                         [0,np.cos(-self.INCLINATION),np.sin(-self.INCLINATION)],
                         [0,-np.sin(-self.INCLINATION),np.cos(-self.INCLINATION)]])
    def R_RA_OF_ASC_NODE(self):
        return np.array([[np.cos(-self.RA_OF_ASC_NODE),np.sin(-self.RA_OF_ASC_NODE),0],
                         [-np.sin(-self.RA_OF_ASC_NODE),np.cos(-self.RA_OF_ASC_NODE),0],
                         [0,0,1]])
    def R_ARG_OF_PERICENTER(self):
        return np.array([[np.cos(-self.ARG_OF_PERICENTER),np.sin(-self.ARG_OF_PERICENTER),0],
                         [-np.sin(-self.ARG_OF_PERICENTER),np.cos(-self.ARG_OF_PERICENTER),0],
                         [0,0,1]])

    def calculate_GMST(self):
        # Calculate the Julian Date (JD) for the given date and time
        a = (14 - self.UPDATE_DATE_MONTH) // 12
        y = self.UPDATE_DATE_YEAR + 4800 - a
        m = self.UPDATE_DATE_MONTH + 12 * a - 3
        JD = (self.UPDATE_DATE_DAY + ((153 * m + 2) // 5) + 365 * y + (y // 4) - (y // 100) + (y // 400) - 32045 +
            (self.UPDATE_DATE_HOUR - 12) / 24.0 + self.UPDATE_DATE_MINUTE / 1440.0 + self.UPDATE_DATE_SECOND / 86400.0)
        # Calculate the Julian centuries since J2000.0
        T = (JD - 2451545.0) / 36525.0
        # Calculate the mean sidereal time in degrees
        GMST = (280.46061837 +
                360.98564736629 * (JD - 2451545.0) +
                T ** 2 * (0.000387933 - T / 38710000))
        # Ensure the result is in the range [0, 360] degrees
        GMST %= 360

        return GMST

  
    def conversion_to_cartesian(self,t):
        self.MEAN_MOTION=self.mean_motion()
        self.MEAN_ANOMALY_UPDATED=self.mean_anomaly(t)
        self.ECCENTRIC_ANOMALY=self.eccentric_anomaly()
        r=self.SEMI_MAJOR_AXIS*(1-self.ECCENTRICITY*np.cos(self.ECCENTRIC_ANOMALY))
        x=self.SEMI_MAJOR_AXIS*(np.cos(self.ECCENTRIC_ANOMALY)-self.ECCENTRICITY)
        y=self.SEMI_MAJOR_AXIS*np.sqrt(1-self.ECCENTRICITY**2)*np.sin(self.ECCENTRIC_ANOMALY)
        x_dot=-np.sin(self.ECCENTRIC_ANOMALY)*self.MEAN_MOTION*self.SEMI_MAJOR_AXIS**2/r
        y_dot=np.sqrt(1-self.ECCENTRICITY**2)*np.cos(self.ECCENTRIC_ANOMALY)*self.MEAN_MOTION*self.SEMI_MAJOR_AXIS**2/r
        pos=self.R_RA_OF_ASC_NODE()@self.R_INCLINATION()@self.R_ARG_OF_PERICENTER()@np.array([[x],[y],[0]])
        vit=self.R_RA_OF_ASC_NODE()@self.R_INCLINATION()@self.R_ARG_OF_PERICENTER()@np.array([[x_dot],[y_dot],[0]])
        fi_earth=t*2*np.pi*(1 + 1/365.25)/(3600*24)
        Rot = np.array([[np.cos(fi_earth), np.sin(fi_earth),0],
                        [-np.sin(fi_earth),np.cos(fi_earth),0],
                        [0,0,1]])
        pos=Rot@pos
        vit=Rot@pos
        return pos.ravel(),vit.ravel()
    
    def conversion_to_geodic_(self,R):
        if len(R.shape)==1:
            R=R.reshape((3,1))
        GMST = self.calculate_GMST()
        p=np.linalg.norm(R,axis=0)
        x=R[0,:]
        y=R[1,:]
        z=R[2,:]
        lon=np.rad2deg(np.arctan2(y,x))-GMST
        lat=np.rad2deg(np.arcsin(z/p))
        return lon, lat
    
    def conversion_to_geodic(self,R,centered=False):
        if len(R.shape)==1:
            R=R.reshape((3,1))
        p=np.linalg.norm(R,axis=0)
        GMST = self.calculate_GMST()
        x = self.R[0,:]
        y = self.R[1,:]
        z = self.R[2,:]
        lon=np.rad2deg(np.arctan2(y,x))-GMST
        lat=np.rad2deg(np.arcsin(z/p))
        for j in range(len(lon)):
            if lon[j] > 180:
                lon[j] -= 360
            elif lon[j] < -180:
                lon[j] += 360
            a = 6378137
            b = 6356752.314235
            f = (a - b) / a

            lat[j] = lat[j] * (1 - f * f)

        return lon, lat
    
    def total(self, nb_orbit=1):
        self.PERIOD=(1/self.MEAN_MOTION)*24*60*60
        self.MEAN_MOTION_SI=2*np.pi/self.PERIOD
        self.SEMI_MAJOR_AXIS=((self.MUE/(self.MEAN_MOTION_SI**2))**(1/3))*1000
        nb_its = nb_orbit*self.PERIOD
        nb_pts = 1000
        T=np.arange(0, nb_its, nb_its/nb_pts)

        self.d1=datetime(self.UPDATE_DATE_YEAR,
                    self.UPDATE_DATE_MONTH,
                    self.UPDATE_DATE_DAY,
                    self.UPDATE_DATE_HOUR,
                    self.UPDATE_DATE_MINUTE,
                    self.UPDATE_DATE_SECOND)
        self.d2=datetime.now(timezone.utc)
        dt=self.d2-self.d1.replace(tzinfo=pytz.UTC)
        
        self.R=np.zeros((3,nb_pts))
        self.Rdot=np.zeros((3,nb_pts))
        for j in range(nb_pts):
            self.R[:,j],self.Rdot[:,j]=self.conversion_to_cartesian(T[j]+float(dt.seconds))
        return self.R
    
    def info(self):
        print("OBJECT_NAME :",self.NAME)
        print("MEAN_MOTION :",self.MEAN_MOTION)
        print("ECCEN :",self.ECCENTRICITY)
        print("INCLI (deg):",np.rad2deg(self.INCLINATION))
        print("RA_NODE :",self.RA_OF_ASC_NODE)
        print("ARG_PERI :",self.ARG_OF_PERICENTER)
        print("MEAN_ANO :",self.MEAN_ANOMALY)
        print("lat :",self.lat[0])
        print("lon :",self.lon[0])

    def plot(self):
        self.lon,self.lat=self.conversion_to_geodic(self.total(),centered=True)
        lon_plot = self.lon.tolist()
        lat_plot = self.lat.tolist()

        for ind_i in range(0,len(lon_plot)-1):
            if np.abs(lon_plot[ind_i]+lon_plot[ind_i+1])<np.abs(lon_plot[ind_i]):
                lon_plot[ind_i]=np.nan
                lat_plot[ind_i]=np.nan
        #self.lon[np.abs(self.lat+90) <= np.abs(self.lat[1]-self.lat[0])] = np.nan
        self.lon,self.lat = lon_plot,lat_plot
        # miller projection 
        map = Basemap(projection='mill',lon_0=0)
        # plot coastlines, draw label meridians and parallels.
        #img = Image.open('map.tif')
        #img_mir = ImageOps.flip(img)
        #map.imshow(img)
        #map.bluemarble()
        
        map.drawcoastlines(color="yellow")
        map.drawparallels(np.arange(-90,90,30),labels=[1,0,0,0])
        map.drawmeridians(np.arange(map.lonmin,map.lonmax+30,60),labels=[0,0,0,1])
        # fill continents 'coral' (with zorder=0), color wet areas 'aqua'
        map.drawmapboundary(fill_color='lightsteelblue')
        map.fillcontinents(color='darkseagreen',lake_color='lightsteelblue')
        map.drawcountries(color="yellow")
        # shade the night areas, with alpha transparency so the 
        # map shows through. Use current time in UTC.
        CS=map.nightshade(self.d2)
        x, y = map(self.lon, self.lat)
        x_0, y_0 = map(self.lon[0], self.lat[0])
        plt.plot(x, y, color="r")
        plt.plot(x_0, y_0, color="r", marker="o")
        plt.title(f"Ground track {self.NAME}")
        plt.show()
        
data = pd.read_json("data.json")
data_sat = data.to_numpy()
data_ID = data["NORAD_CAT_ID"].to_numpy()
input_name = int(input("Enter the NORAD ID: "))#
index =  np.where(data_ID == input_name)[0][0]
print(f"You chose: {data_sat[index][0]}")
satellite = Satellite_class(data_sat[index])
satellite.plot()


#satellite.info()
