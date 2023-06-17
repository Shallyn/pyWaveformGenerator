"""
This is the module for gravitational wave coherent search.
Writer: Shallyn(shallyn.liu@foxmail.com)
"""

import numpy as np
import astropy.time as at
from astropy.constants import c as astro_c
import ephem
c_SI = astro_c.value # m s-1
#-------------------GMST-------------------#
def gmst_accurate(gps_time):
    # gmst = at.Time(gps_time, format='gps',
    #             location=(0, 0)).sidereal_time('mean').rad
    obs = ephem.Observer()
    obs.lon, obs.lat = 0, 0
    gps = at.Time(gps_time, format = 'gps')
    obs.date = gps.datetime
    ret = 1.*obs.sidereal_time()
    return ret


def time_delay(ifo, ra, de, gpstime):
    det = Detector(ifo)
    return det.time_delay_from_earth_center(ra, de, gpstime)

class ifoinfo(object):
    def __init__(self, name, location, latitude, longtitude):
        self.location = np.array(location)
        self.latitude = latitude
        self.longtitude = longtitude
    
    def setARMResponse(self, xdx, xdy, xdz, ydx, ydy, ydz):
        response = np.zeros([3,3])
        response[0,0] = (xdx **2 - ydx **2)/2
        response[0,1] = (xdx*xdy - ydx*ydy)/2
        response[0,2] = (xdx*xdz - ydx*ydz)/2
        
        response[1,0] = (xdy*xdx - ydy*ydx)/2
        response[1,1] = (xdy **2 - ydy **2)/2
        response[1,2] = (xdy*xdz - ydy*ydz)/2
        
        response[2,0] = (xdz*xdx - ydz*ydx)/2
        response[2,1] = (xdz*xdy - ydz*ydy)/2
        response[2,2] = (xdz **2 - ydz **2)/2
        self.response = response

# distance in [m]
#H1
H1 = ifoinfo('H1', \
             location = [-2.16141492636e+06, \
                         -3.83469517889e+06, \
                         4.60035022664e+06], \
             latitude = 0.81079526383, \
             longtitude = -2.08405676917)
H1.setARMResponse(xdx = -0.22389266154,\
                  xdy = 0.79983062746,\
                  xdz = 0.55690487831,\
                  ydx = -0.91397818574,\
                  ydy = 0.02609403989,\
                  ydz = -0.40492342125)
#L1
L1 = ifoinfo('L1', \
             location = [-7.42760447238e+04, \
                         -5.49628371971e+06, \
                         3.22425701744e+06], \
             latitude = 0.53342313506, \
             longtitude = -1.58430937078)
L1.setARMResponse(xdx = -0.95457412153,\
                  xdy = -0.14158077340,\
                  xdz = -0.26218911324,\
                  ydx = 0.29774156894,\
                  ydy = -0.48791033647,\
                  ydz = -0.82054461286)
#V1
V1 = ifoinfo('V1', \
             location = [4.54637409900e+06, \
                         8.42989697626e+05, \
                         4.37857696241e+06], \
             latitude = 0.76151183984, \
             longtitude = 0.18333805213)
V1.setARMResponse(xdx = -0.70045821479,\
                  xdy = 0.20848948619,\
                  xdz = 0.68256166277,\
                  ydx = -0.05379255368,\
                  ydy = -0.96908180549,\
                  ydz = 0.24080451708)

#K1
K1 = ifoinfo('K1', \
             location = [-3777336.024, \
                         3484898.411, \
                         3765313.697], \
             latitude = 0.6355068497, \
             longtitude = 2.396441015)
K1.setARMResponse(xdx = -0.3759040,\
                  xdy = -0.8361583,\
                  xdz = 0.3994189,\
                  ydx = 0.7164378,\
                  ydy = 0.01114076,\
                  ydz = 0.6975620)

#ET1
E1 = ifoinfo('E1', \
             location = [4.54637409900e+06, \
                         8.42989697626e+05, \
                         4.37857696241e+06], \
             latitude = 0.76151183984, \
             longtitude = 0.18333805213)
E1.setARMResponse(xdx = -0.70045821479,\
                  xdy = 0.20848948619,\
                  xdz = 0.68256166277,\
                  ydx = -0.39681482542,\
                  ydy = -0.73500471881,\
                  ydz = 0.54982366052)

#ET2
E2 = ifoinfo('E2', \
             location = [4.53936951685e+06, \
                         8.45074592488e+05, \
                         4.38540257904e+06], \
             latitude = 0.76299307990, \
             longtitude = 0.18405858870)
E2.setARMResponse(xdx = 0.30364338937,\
                  xdy = -0.94349420500,\
                  xdz = -0.13273800225,\
                  ydx = 0.70045821479 ,\
                  ydy = -0.20848948619,\
                  ydz = -0.68256166277)

#ET3
E3 = ifoinfo('E3', \
             location = [4.54240595075e+06, \
                         8.35639650438e+05, \
                         4.38407519902e+06], \
             latitude = 0.76270463257, \
             longtitude = 0.18192996730)
E3.setARMResponse(xdx = 0.39681482542,\
                  xdy = 0.73500471881,\
                  xdz = -0.54982366052,\
                  ydx = -0.30364338937 ,\
                  ydy = 0.94349420500,\
                  ydz = 0.13273800225)

_ifodict = dict()
_ifodict['H1'] = H1
_ifodict['L1'] = L1
_ifodict['V1'] = V1
_ifodict['K1'] = K1
_ifodict['E1'] = E1
_ifodict['E2'] = E2
_ifodict['E3'] = E3

class Detector(object):
    def __init__(self, name):
        if name not in _ifodict:
            msg = 'No such ifo: {}'.format(name)
            raise ValueError(msg)
        self.frDetector = _ifodict[name]
        self.response = self.frDetector.response
        self.location = self.frDetector.location
        self.latitude = self.frDetector.latitude
        self.longtitude = self.frDetector.longtitude

    def get_at_and_delay(self, ra, de, psi, gps):
        gmst = gmst_accurate(gps)
        delay = self.time_delay_from_earth_center_gmst(ra, de, gmst)
        at = self.antenna_pattern_gmst(ra, de, psi, gmst)
        return (at, delay)
    
    def time_delay_from_earth_center(self, ra, de, gps):
        gcloc = np.array([0,0,0])
        gha = gmst_accurate(gps) - ra
        ehat_src = np.zeros(3)
        ehat_src[0] = np.cos(de) * np.cos(gha)
        ehat_src[1] = np.cos(de) * (-np.sin(gha))
        ehat_src[2] = np.sin(de)
        delta_xyz = gcloc - self.location
        delta_xyz.reshape(gcloc.shape)
        return np.dot(ehat_src, delta_xyz) / c_SI
    
    def time_delay_from_earth_center_gmst(self, ra, de, gmst):
        gcloc = np.array([0,0,0])
        gha = gmst - ra
        ehat_src = np.zeros(3)
        ehat_src[0] = np.cos(de) * np.cos(gha)
        ehat_src[1] = np.cos(de) * (-np.sin(gha))
        ehat_src[2] = np.sin(de)
        delta_xyz = gcloc - self.location
        delta_xyz.reshape(gcloc.shape)
        return np.dot(ehat_src, delta_xyz) / c_SI

    
    def antenna_pattern(self, ra, de, psi, gps):
        gmst = gmst_accurate(gps)
        D = self.response
        gha = gmst - ra
        cosgha = np.cos(gha)
        singha = np.sin(gha)
        cosdec = np.cos(de)
        sindec = np.sin(de)
        cospsi = np.cos(psi)
        sinpsi = np.sin(psi)

        x0 = -cospsi * singha - sinpsi * cosgha * sindec
        x1 = -cospsi * cosgha + sinpsi * singha * sindec
        x2 =  sinpsi * cosdec
        x = np.array([x0, x1, x2])

        dx = np.dot(D, x)

        y0 =  sinpsi * singha - cospsi * cosgha * sindec
        y1 =  sinpsi * cosgha + cospsi * singha * sindec
        y2 =  cospsi * cosdec
        y = np.array([y0, y1, y2])
        dy = np.dot(D, y)
        Fplus = (x * dx - y * dy).sum(axis=0)
        Fcross = (x * dy + y * dx).sum(axis=0)
        return Fplus, Fcross

    def antenna_pattern_gmst(self, ra, de, psi, gmst):
        D = self.response
        gha = gmst - ra
        cosgha = np.cos(gha)
        singha = np.sin(gha)
        cosdec = np.cos(de)
        sindec = np.sin(de)
        cospsi = np.cos(psi)
        sinpsi = np.sin(psi)

        x0 = -cospsi * singha - sinpsi * cosgha * sindec
        x1 = -cospsi * cosgha + sinpsi * singha * sindec
        x2 =  sinpsi * cosdec
        x = np.array([x0, x1, x2])

        dx = np.dot(D, x)

        y0 =  sinpsi * singha - cospsi * cosgha * sindec
        y1 =  sinpsi * cosgha + cospsi * singha * sindec
        y2 =  cospsi * cosdec
        y = np.array([y0, y1, y2])
        dy = np.dot(D, y)
        Fplus = (x * dx - y * dy).sum()
        Fcross = (x * dy + y * dx).sum()
        return Fplus, Fcross
    
    def amplitude_modulation(self, ra, de, gmst):
        gha = gmst - ra
        D = self.response
        cosgha = np.cos(gha)
        singha = np.sin(gha)
        cosdec = np.cos(de)
        sindec = np.sin(de)

        x0 = -singha
        x1 = -cosgha
        x2 = 0
        x = np.array([x0, x1, x2])

        dx = np.dot(D, x)

        y0 =  -cosgha * sindec
        y1 =  singha * sindec
        y2 =  cosdec
        y = np.array([y0, y1, y2])
        dy = np.dot(D, y)
        Gplus = (x * dx - y * dy).sum()
        Gcross = (x * dy + y * dx).sum()
        return Gplus, Gcross
