import numpy as np
from astropy.coordinates import *
from astropy.time import Time
import astropy.units as u

time = Time("2024-10-23 17:00:00")
sun = get_sun(time)

observer =  EarthLocation(lat=38.8976841639134 * u.deg, lon=-77.03653291227445 * u.deg)
sunaltaz = sun.transform_to(AltAz(obstime=time, location=observer))

phi_sun = np.deg2rad(sunaltaz.az)
altitude_sun = np.deg2rad(sunaltaz.alt)
theta_sun = np.deg2rad(90*u.deg) - altitude_sun